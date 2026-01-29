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
    Foam::constantLhgr

Description
	The `constantLhgr` heatSource class in OFFBEAT allows you to set a constant
	average **linear heat generation rate (lhgr)** via a keyword in the heatSource
	subDictionary.

	It is also possible to apply a radial and an axial profile to the lhgr. Note
	that even if the lhgr is constant over time, the axial and radial profile might
	change over time, depending on the profile type chosen by the user.
	
	#### Formulation
	
	The final volumetric heat source is derived as a combination of lhgr, radial and
	axial profiles, as follows:

	$$
	\begin{aligned}
	Q(t,r,z) = Q_{avg} \cdot f(r,t) \cdot g(z,t)
	\end{aligned}
	$$
	
	where:
	
	- \f$r\f$ is the relative radial position (0...1)
	- \f$z\f$ is the relative axial position (0...1)
	- \f$t\f$ is the current time
	- \f$f(r,t)\f$ and \f$g(z,t)\f$ are the radial and axial profiles, respectively
	- \f$Q_{avg}\f$ is the average power density in \f$W/m^3\f$.
	
	The radial and axial profiles are expected to be normalized to 1, while the
	average power density is obtained from:

	$$
	\begin{aligned}
	Q_{avg} = \frac{lhgr \cdot h_{ref}}{\frac{V_{model}}{\alpha}}
	\end{aligned}
	$$
	
	where:
	
	- \f$h_{ref}\f$ is the height of the rod in the reference, undeformed geometry
	- \f$V_{model}\f$ is the volume of the fuel as calculated by OpenFOAM in the
	model
	- \f$\alpha\f$ is the fraction of the real geometry modeled in the simulations
	(e.g. 0.25 if you are using symmetry BCs and you are simulating only a quarter
	of the rod).
	
________________________________________________________________________________

\par Options
	The `constantLhgr` heat source model requires the user to specify a few
	additional parameters in the `heatSourceOptions` sub-dictionary located 
	within the `solverDict`.
    
    Parameters in `heatSourceOptions`:
	
	- <b>`lhgr`</b> - The constant linear heat generation rate in W/m.
	- <b>`materials`</b> - A list of material (or cellZones names) where the heat source
	model applies.
	- <b>`axialProfile`</b> - A sub-dictionary that specifies the type of axial profile.
	- <b>`radialProfile`</b> - A sub-dictionary that specifies the type of radial profile.
	- <b>`useExternalReferenceDimensions`</b> – A bool to activate reading reference dimensions for 
	axial profile (`zMin` and `zMax`) from the dictionary (false by default).  
	- <b>`zMin`</b> – The reference axial coordinate of the bottom of the fuel region. 
	Only used if `useReferenceDimensions` is set to true.  
	- <b>`zMax`</b> – The reference axial coordinate of the top of the fuel region. 
	Only used if `useReferenceDimensions` is set to true.
	- <b>`volumeFraction`</b> – The fraction of the reference rod present in the local 
	mesh. Only used if `useReferenceDimensions` is set to true.

	
	The `constantLhgr` heat source model is designed to calculate the volumetric
	power density for a cylindrical rod with an initial reference height (i.e., in
	undeformed geometry) along a given axial direction. For simulations where there
	is an axis of symmetry, the class can take into consideration the angular
	fraction of the real rod that is being modeled.
	
	Thus, the class requires the additional parameters in the `globalOptions` sub-dictionary 
	of the `solverDict` dictionary.
	
    Parameters in `globalOptions`:

	- <b>`pinDirection`</b> - The vector of the axial direction (e.g. (0 0 1) for the
	z-axis).
	- <b>`angularFraction`</b> - The fraction of the 360 degree angle that is being
	modeled. 
	
	\warning
	Reove the <b>`angularFraction`</b> keyword for axisymmetric r-z models, as the 
	angular fraction is automatically calculated from the wedge boundaries. Otherwise, it
	must be provided. For example, set it to 0.25 for quarter symmetry, to 0.5 for
	half-symmetry and to 1.0 for 3D simulations.
	
________________________________________________________________________________

Usage
	Here is a code snippet the `solverDict` to be used for the case of a constant
	lhgr of 10e3 W/m with flat (i.e. uniform) axial and radial profiles.
	This particular example considers a quarter of a fuel disc.
	
	\verbatim
		heatSource constantLhgr;

		globalOptions
		{
				angularFraction 0.25;
				pinDirection  (0 0 1);
		}

		heatSourceOptions
		{
			lhgr        10e3;

			materials ( fuel );

			useExternalReferenceDimensions off;
			// // Alternatively one can use external reference dimensions
			// useExternalReferenceDimensions on;
			// zMin 0.0;
			// zmax 3.0;
			// volumeFraction 0.01;

			axialProfile
			{
				type flat;
			}

			radialProfile
			{
				type flat;
			}
		}
	\endverbatim
	
	For a time-dependent axial profile, one can modify the `axialProfile` sub-
	dictionary as follows:
	
	\verbatim
	heatSourceOptions
	{
		// .. The rest of the dictionary may remain unchanged

		axialProfile
		{
			// Example: time dependent axial radial profile.
			// The user provides a 2D table (z-vs-time) of normalized lhgr values
			// (normalized to the average).
			// Other options are available (see axialProfile class).
			type            timeDependentTabulated;

			// List of time points for axial profile
			timePoints      ( 0 3600 31539600 );

			// List of axial locations for axial profile
			axialLocations  ( 0.0 0.5 1.0 );

			// 2D table for axial profile (z is x axis; time is y axis)
			data
			(
				( 1.0 1.0 1.0 )
				( 0.5 1.0 1.5 )
				( 0.8 1.1 1.0 )
			);

			// Select the axial interpolation method
			axialInterpolationMethod    linear;

			// Select the time interpolation method for time steps that fall in-between
			// the time points indicated in the timePoints list above
			timeInterpolationMethod     linear;
		}

	}
	\endverbatim

SourceFiles
    heatSource.C

\todo
    - The lhgr might be needed by some models such as the relocation, therefore
    it might be useful to have it as a field in the registry. Alternatively, 
    one could make the heatSource class an IOObject and ask directly for the lhgr 
    with a function. In this case having the lhgr as a vol field might not be necessary.

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef constantLhgr_H
#define constantLhgr_H

#include "constantHeatSource.H"
#include "axialProfile.H"
#include "radialProfile.H"
#include "azimuthalProfile.H"
#include "volFields.H"
#include "InterpolateTables.H"
#include "globalOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class constantLhgr Declaration
\*---------------------------------------------------------------------------*/

class constantLhgr
:
    public constantHeatSource
{
 
    // Private data
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        constantLhgr(const constantLhgr&);

        //- Disallow default bitwise assignment
        void operator=(const constantLhgr&);
        
protected:   

    // Protected data
        
        //- Linear heat generation rate (W/m)
        volScalarField lhgr_;

        //- Constant lhgr average value
        scalar constLhgr_;
        
        //- Pin axial direction
        vector pinDirection_;
        
        //- Radial direction
        vector radialDirection_;
        
        //- Minimum z
        scalar zMin_;
        
        //- Maximum z
        scalar zMax_;

		//- Bool to read reference dimensions from ditionary
		Switch useExternalReferenceDimensions_;

		//- Fraction of the reference volume (e.g., full rod) present in the mesh
		scalar refVolumeFraction_;

        //- Maximum radius
        scalar rMax_;

        //- Time-dependent radial profile
        autoPtr<radialProfile> radialModel_;

        //- Time-dependent azimuthal profile
        autoPtr<azimuthalProfile> azimuthalModel_;
        
        //- Time-dependent axial profile
        autoPtr<axialProfile> axialModel_;

        //- Addressing for power-producing cells
        labelList addr_;

public:

    //- Runtime type information
    TypeName("constantLhgr");

    // Constructors

        //- Construct from mesh, materials and dict
        constantLhgr
        (
            const fvMesh& mesh,
            const materials& materials,
            const dictionary& heatSourceDict,
            const bool timeDependent=false
        );


    //- Destructor
    ~constantLhgr();


    // Member Functions

        //- Get addressing for power-producing cells
        void calcAddressing(const dictionary& dict);
        
        //- Calculate the fuel height
        void calcReferenceDimensions(const dictionary& dict);
    
        //- update the power distribution in the supplied list of cells
        virtual void correct();
        
        //- Return the next average power density (assuming the deltaT stays the
        // same)
        virtual scalar predictNextPowerDensity() const;    

        //- Return the average power density 
        virtual scalar averagePowerDensity() const;

        //- Reshape Q profile according to porosity
        virtual void reshapePowerPorosity();

        //- For time-dependent heat sources, return the next time marker,
        //- i.e. the next time point at which the rate of power change
        //- will be discontinuous or non-smooth.
        virtual scalar nextTimeMarker() const;

        //- Return the final time marker after which the power distribution
        //- is no longer defined.
        virtual scalar lastTimeMarker() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
