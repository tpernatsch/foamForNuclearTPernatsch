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
    Foam::burnupLassmannFBR

Description
	The `LassmannFBR` burnup class in OFFBEAT is inspired by the work of Lassmann et
	al., "The radial distribution of plutonium in high burnup UO2 fuels" Journal of
	Nuclear Materials.

	The class extends the functionality of the standard `Lassmann` burnup class 
    to fast spectrum applications. In addition to the standard burnup class, 
    this extension allows for:

    
    - Solving for depletion of an extended number of nuclides (He4, U234, U235, 
    U236, U237, U238, Np237, Np238, Np239, Pu238, Pu239, Pu240, Pu241, Pu242, 
    Pu243, Am241, Am242, Am242m, Am243, Am244, Cm242, Cm243, Cm244 and Cm245);

    - Accounting for He production through the mechanisms of alpha-decay, 
    ternary fission and (n/alpha) reaction of O-16;
    
    - Reading cross section data from lookup tables located in `constant/XSdata`;

    The `constant/XSdata` dictionary allows the user to provide parametrized 
    cross sectionf to be used for nuclide depletion. The type of parametrization
    depends on the `reactorType` selected by the user:

    - If `LWR` is selected, fission and capture cross section are parametrized 
    with respect the U-235 enrichment in UO2;

    - If `FBR` is selected, fission and capture cross section are parametrized 
    with respect the Pu wt% concentration in MOX;

    \warning
	The models requires the user to specify the `reactorType` keyword in the 
    `globalOptions` dictionary of `coinstant/solverDict`. Currently, only `LWR` 
    and `FBR` reactor types are supported.
	
________________________________________________________________________________

\par Options	

    Parameters in `materials` (for all fuel-type materials):
	
	- <b>`enrichment`</b> - For UO2 fuels, the enrichment of fissile material 
    in fraction (i.e. 0.025 for a 2.5% enrichment).
	- <b>`isotopes`</b> - For MOX fuels, the isotopic composition of the fuel
	- <b>`densityFraction`</b> - The fraction of theoretical density of the fuel (e.g.
	0.95).
	
	Finally, the `LassmannFBR` burnup model requires the use of a `neutronicsSolver`
	(for instance, the `diffusion` solver) that calculates the neutron flux
	distribution in the fuel.
	
________________________________________________________________________________

Usage
	Here is an example of the `solverDict` to be used for activating the 
    `LassmannFBR` burnup model
	
	\verbatim
        globalOptions
        {
            pinDirection            (0 0 1);
            reactorType             "FBR";
        }

        burnup LassmannFBR;
        neutronicsSolver diffusion;

        // ... 
	
        materials
        {
            fuel
            {
                material UPuO2;

                oxygenMetalRatio            1.97;

                isotopes
                {
                    U
                    {
                        ratioOverMetal 0.80;
                        weightFractions
                        (
                            235 0.007
                            238 0.993
                        );
                    }
                    Pu
                    {
                        ratioOverMetal 0.20;
                        weightFractions 
                        ( 
                            238 2.332e-2 
                            239 56.873e-2 
                            240 26.997e-2 
                            241 6.105e-2 
                            242 7.693e-2 
                        );
                    }
                }

            }

                // ... Rest of `fuel` dictionary
            }

            //.. Rest of materials
        
        }
	\endverbatim

SourceFiles
    burnupLassmannFBR.C

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

#ifndef burnupLassmannFBR_H
#define burnupLassmannFBR_H

#include "volFields.H"
#include "objectRegistry.H"
#include "IOdictionary.H"

#include "fuelMaterial.H"
#include "UO2.H"
#include "UPuO2.H"
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
                         Class burnupLassmannFBR Declaration
\*---------------------------------------------------------------------------*/

class burnupLassmannFBR
:
    public burnupFromPower
{
    // Private data

        //- Const ref to global options
        const globalOptions& globOpt_;

        //- Dictionary for XS 
        IOdictionary XSdict_;

        //- Power density form factor
        volScalarField formFactor_;

        //- Burnup in MWd/kgU
        volScalarField BukgU_;

        //- HashTable with unique label for each nuclide
        HashTable<label> ID_;

        //- HashTable for nuclides decay coefficient
        HashTable<scalar> nuclLambdas_;

        //- Ptr list of volScalarFields for nuclides conc. [1/cm3]
        PtrList<volScalarField> nuclPtrList_;

        //- Ptr list of volScalarFields for elements conc (tot. U, Pu) [1/cm3]
        PtrList<volScalarField> elemPtrList_;

        //- Microscopic fission cross section
        scalarField sf_;

        //- Microscopic capture cross section
        scalarField sc_;

        //- Microscopic absorption cross section
        scalarField sa_;

        //- Thermal absorption microscopic cross section
        scalarField sa_th_;

        //- XS for U238(n,2n) reaction
        scalar sn2nU238_;

        //- XS for U238(n,3n) reaction
        scalar sn3nU238_;

        //- XS for O16 (n,alpha) reaction
        scalar snAlphaO16_;

        //- Branching ratio beta decay Am242
        scalar BR_beta_Am242_;

        //- Branching ratio for Am-241(n,g)Am-242
        scalar BR_ng_Am242_;

        //- Isomeric Branching Ratio Am241 (g)
        scalar BR_Am241_;

        //- Isomeric Branching Ratio Am243 (g)
        scalar BR_Am243_;

        //- Ratio of Heavy Metals over fuel (assumed 0.8815)
        const scalar HMoverFuelRatio_;

        //- Mean energy per fission event in (MWd)
        scalar Efiss_;

        //- Pellet inner radius (for each cell)
        scalarField rIn_;

        //- Pellet outer radius (for each cell) 
        scalarField rOut_;

        //- Convergence criterion
        scalar precision_;

        //- Mapper for averaging quantities in axial slices
        const sliceMapper* mapper_;

    // Private Member Functions

        //- Set hash table of nuclides IDs
        virtual HashTable<label> setNuclidesIDs(); 

        //- Set hash table of nuclides decay constants
        virtual HashTable<scalar> setNuclidesLambdas(); 

        // Fill pointer list with nuclides volScalarFields
        virtual void createNuclidesFields();

        //- Set cross sections of each nuclides for each material
        virtual void setXS(); 

        //- Disallow default bitwise copy construct
        burnupLassmannFBR(const burnupLassmannFBR&);

        //- Disallow default bitwise assignment
        void operator=(const burnupLassmannFBR&);

public:

    //- Runtime type information
    TypeName("LassmannFBR");

    // Constructors

    //- Construct from mesh, global Options, materials and dict
    burnupLassmannFBR
    (
        const fvMesh& mesh,
        const materials& mat,
        const dictionary& burnupDict
    );

    //- Destructor
    ~burnupLassmannFBR();

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

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
