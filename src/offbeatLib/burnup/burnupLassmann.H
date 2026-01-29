/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::burnupLassmann

Description
	The `Lassmann` burnup class in OFFBEAT is inspired by the work of Lassmann et
	al., "The radial distribution of plutonium in high burnup UO2 fuels" Journal of
	Nuclear Materials.

	The class extends the functionality of the base `fromPower` burnup class as it
	serves as a simplified depletion module incorporating a small set of relevant
	nuclides and solving the Bateman equation for their time evolution.
	
	The class is designed to provide a more accurate representation of the power
	distribution within the fuel by considering the spatial variation of power
	density across the radial domain.
	
	\note
	The `Lassmann` burnup model is designed to be used with axially symmetric models
	or 3D models where there is no significant change in power/nuclide along the
	azimuthal angle. For each slice, the form factor in the fuel will only be 
    calculated along the radius, with no variation alon the azimuthal angle.

    \warning
	The models require the subdivision of the fuel rod into axial slices. Thus
	one should activate the `sliceMapper` model in OFFBEAT. Also, the `neutronicsSolver` 
    should be activated as a way of providing the neutron flux distribution to 
    the depletion module.
	
	#### Formulation
	
	First, at each time step, the updated average burnup for each fuel slice is
	calculated as:
	
	$$
	\begin{aligned}
	Bu = Bu^0 + \frac{Q \cdot \Delta t}{1000 \cdot \rho}
	\end{aligned}
	$$
	
	where:
	
	- \f$Bu^0\f$ is the old local burnup value, i.e. at the end of the previous time
	step, in MWd/MT\f$_{oxide}\f$
	- \f$Q\f$ is the power density in W/m\f$^3\f$
	- \f$\rho\f$ is density in kg/m\f$^3\f$
	- \f$\Delta t\f$ is time step in seconds
	
	Then, a simplified version of the Bateman equations is solved for a small set of
	relevant nuclides. At the moment: U235, U236, U238, Np237, Pu238, Pu239, Pu240,
	Pu241, Pu242, Am241, Am243, Cm242, Cm243, Cm244.
	These nuclides are chosen based on their significance for the power profile.
	
	The updated values of the fissile nuclide distribution in each cell are combined
	with the flux distribution obtained from the `neutronicsSolver` model to obtain
	the power profile in a given slice:
	
	$$
	\begin{aligned}
	Q(r) = \sum_f \left( {N_f(r) \sigma_f(r)} \phi(r) \right)
	\end{aligned}
	$$
	
	where \f$\phi\f$ is the flux, N\f$_f\f$ is the nuclide concentration of a
	fissile isotope and \f$\sigma_f\f$ is the fissile cross section,
	
	The local power density is then used to calculate the local burnup increment,
	which quantifies the amount of fuel burnup that occurs in a given region during
	a specific time step.

	This increment is used to update the nuclide distribution for the next iteration
	of the Lassman model.
	
	For a detailed explanation of the Bateman equation and its application in burnup
	calculations, it is recommended to refer to relevant nuclear engineering
	literature or online resources.

________________________________________________________________________________

\par Options	
	The `Lassmann` burnup model has at the moment only a single dedicated parameter
	that the user can specify in the `burnupOptions` sub-dictionary located within the
    `solverDict`. In addition, the `Lassmann` burnup model requires that the materials	
    sub-dictionaries (always within the `solverDict`) for all fuel-type materials such 
    as U2 or MOX, include some additional parameters.
	
    Parameters in `burnupOptions`:

    - <b>`convergePrecision`</b> - Parameter that controls the convergence precision for
	the burnup calculations. **By default, it is set to 1e-2 (which is the suggested
	value).**
    
    Parameters in `materials` (for all fuel-type materials):
	
	- <b>`enrichment`</b> - The enrichment of fissile material, in fraction (i.e. 0.025 for
	a 2.5% enrichment).
	- <b>`densityFraction`</b> - The fraction of theoretical density of the fuel (e.g.
	0.95).
	
	Finally, the `Lassmann` burnup model requires the use of a `neutronicsSolver`
	(for instance, the `diffusion` solver) that calculates the neutron flux
	distribution in the fuel.
	
________________________________________________________________________________

Usage
	Here is an example of the `solverDict` to be used for activating the `Lassmann`
	burnup model
	
	\verbatim
        burnup Lassmann;
        neutronicsSolver diffusion;
        
        burnupOptions
        {
        	// Name of the heat source field (set to "Q" by default)
            heatSourceName  "Q";
        
            // Convergence criterion for depletion iterations. By default set to 1e-2.
            convergencePrecistion 1e-2;
        }

        // ... Rest of solverDict
	
        materials
        {
            fuel
            {
                material UO2;

                enrichment 0.025;
                densityFraction 0.95;

                // ... Rest of `fuel` dictionary
            }

            //.. Rest of materials
        
        }
	\endverbatim

SourceFiles
    burnupLassmann.C

\mainauthor
    Y. Robert - INSA, Lyon\n
    A. Scolaro, E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE,
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef burnupLassmann_H
#define burnupLassmann_H

#include "volFields.H"
#include "objectRegistry.H"
#include "IOdictionary.H"

#include "fuelMaterial.H"
#include "UO2.H"
#include "sliceMapper.H"

#include "primitiveFields.H"
#include "fvMesh.H"
#include "surfaceFields.H"

#include "burnupFromPower.H"
#include "globalOptions.H"

#include "InterpolateTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class burnupLassmann Declaration
\*---------------------------------------------------------------------------*/

class burnupLassmann
:
    public burnupFromPower
{
    // Private data

        //- Power density form factor
        volScalarField formFactor_;

        //- Burnup in MWd/kgU
        volScalarField BukgU_;

        //- Number of elements for depletion
        const scalar Nelements_;

        //- Atomic densities for U235
        volScalarField& N_U235_;

        //- Atomic densities for U236
        volScalarField& N_U236_;

        //- Atomic densities for Np237
        volScalarField& N_Np237_;

        //- Atomic densities for Pu238
        volScalarField& N_Pu238_;

        //- Atomic densities for U238
        volScalarField& N_U238_;

        //- Atomic densities for Pu239
        volScalarField& N_Pu239_;

        //- Atomic densities for Pu240
        volScalarField& N_Pu240_;

        //- Atomic densities for Pu241
        volScalarField& N_Pu241_;

        //- Atomic densities for Pu242
        volScalarField& N_Pu242_;

        //- Atomic densities for Am241
        volScalarField& N_Am241_;

        //- Atomic densities for Am243
        volScalarField& N_Am243_;

        //- Atomic densities for Cm242
        volScalarField& N_Cm242_;

        //- Atomic densities for Cm243
        volScalarField& N_Cm243_;

        //- Atomic densities for Cm244
        volScalarField& N_Cm244_;

        //- Nuclides index for U235
        const int iU235_  ;

        //- Nuclides index for U236
        const int iU236_  ;

        //- Nuclides index for Np237
        const int iNp237_ ;

        //- Nuclides index for Pu238
        const int iPu238_ ;

        //- Nuclides index for U238
        const int iU238_  ;

        //- Nuclides index for Pu239
        const int iPu239_ ;

        //- Nuclides index for Pu240
        const int iPu240_ ;

        //- Nuclides index for Pu241
        const int iPu241_ ;

        //- Nuclides index for Pu242
        const int iPu242_ ;

        //- Nuclides index for Am241
        const int iAm241_ ;

        //- Nuclides index for Am243
        const int iAm243_ ;

        //- Nuclides index for Cm242
        const int iCm242_ ;

        //- Nuclides index for Cm243
        const int iCm243_ ;

        //- Nuclides index for Cm244
        const int iCm244_ ;

        //- Microscopic fission cross section
        scalarField sf_;

        //- Microscopic capture cross section
        scalarField sc_;

        //- Thermal absorption microscopic cross section
        scalarField sa_th_;

        //- XS for U238(n,2n)Np237 reaction
        scalar sn2n_;

        //- Pellet inner radius (for each cell)
        scalarField rIn_;

        //- Pellet outer radius (for each cell) 
        scalarField rOut_;

        //- U235 enrichement
        scalarField enrichment_;

        //- Fuel density fraction    
        scalarField densityFraction_;    

        //- Convergence criterion
        scalar precision_; 

        //- Mapper for averaging quantities in axial slices
        const sliceMapper* mapper_;


    // Private Member Functions
    
        //- Set one-group fission and capture XS in the given fields
        //  for the specified U-235 enrichment (integer wt-%)
        void setXSForIntEnrichment(const int, scalarField&, scalarField&) const;

        //- Disallow default bitwise copy construct
        burnupLassmann(const burnupLassmann&);

        //- Disallow default bitwise assignment
        void operator=(const burnupLassmann&);

public:

    //- Runtime type information
    TypeName("Lassmann");

    // Constructors

    //- Construct from mesh, global Options, materials and dict
    burnupLassmann
    (
        const fvMesh& mesh,
        const materials& mat,
        const dictionary& burnupDict
    );

    //- Destructor
    ~burnupLassmann();

    // Access Functions

    // Member Functions
    
    //- Update the burnup
    virtual void correct();

    //- Return a const reference to SigmaA_ 
    virtual const scalarField SigmaA() const;

    //- Return a const reference to the diffusion coeff field 
    virtual const scalarField neutronD() const;

    //- Calculate inner and outer radii
    virtual void calcRinRout(const labelList& addr);

    //- Initialize XS and nuclear data
    virtual void initNuclearData(const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
