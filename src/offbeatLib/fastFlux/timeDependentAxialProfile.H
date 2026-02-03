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
    Foam::timeDependentAxialProfile

Description
	The `timeDependentAxialProfile` class in OFFBEAT allows you to set a time
	dependent average fast flux.
	
	It is also possible to apply a an axial profile to the fast flux. Note that even
	if the fast flux is constant over time, the axial profile might change over
	time, depending on the profile type chosen by the user.
	
	The fast fluence will evolve in time according to the local value of the fast
	flux.
	
	#### Formulation
	
	The final local fast flux value is derived as a combination of rod-average fast
	flux and axial profile, as follows:
	
	$$
	\begin{aligned}
	fastFlux(t, z) = \phi(t) \cdot g(z,t)
	\end{aligned}
	$$
	
	where:
	
	- \f$z\f$ is the relative axial position (0...1)
	- \f$t\f$ is the current time
	- \f$g(z,t)\f$ is the axial profile
	- \f$\phi\f$ is the rod-average fast flux in <code>n/cm<sup>2</sup>/s</code>.
	
	The axial profile is expected to be normalized to 1.
	
________________________________________________________________________________

\par Options	
	The `timeDependentAxialProfile` fast flux model requires the user to specify a
	few additional parameters in the `fastFluxOptions` sub-dictionary located
	within the `solverDict`.	
    
    Parameters in `fastFluxOptions`:
	
	- <b>`timePoints`</b> - A list of time values at which the \f$\phi\f$ is provided. The time
	unit depend on the userTime selected by the user (seconds by default).
	- <b>`fastFlux`</b> - A list of \f$\phi\f$ values in \f$n/cm^2/s\f$, one value per
	time-point.
	- <b>`timeInterpolationMethod`</b> - Select the time interpolation method for time
	steps that fall in between the time points indicated in the timePoints list.
	- <b>`materials`</b> - A list of material (or cellZones names) where the fast flux
	source model applies.
	- <b>`axialProfile`</b> - A sub-dictionary that specifies the type of axial profile.
	
________________________________________________________________________________

Usage
	Here is a code snippet of the `solverDict` to be used for the case of a fast
	flux that ramps from 0 to 1e13 <code>n/cm<sup>2</sup>/s</code> in 1hr time, remains
	constant for a year, and then ramps down to 0 in 1 hr.
	This particular example considers  flat (i.e. uniform) axial profile.
	
	\verbatim
		fastFlux timeDependentAxialProfile;
		
		fastFluxOptions
		{
			//           t0  1hr    1yr       1yr+1hr
			timePoints  ( 0   3600   31536000  31539600);
			fastFlux    ( 0   1e13   1e13      0       );
			
			timeInterpolationMethod linear; //step;
			
			materials ( fuel cladding);
			
			axialProfile
			{
				type flat;
			}
		}
	\endverbatim
	
	For a time-dependent axial profile, one can modify the `axialProfile` sub-
	dictionary selecting the `timeDependentTabulated`. Other options might be available 
	(see axialProfile class):
	
	\verbatim
		fastFluxOptions
		{
			// .. The rest of the dictionary may remain unchanged
			
			axialProfile
			{
				// Example: time dependent axial radial profile.
				type            timeDependentTabulated;
				
				// List of time points for axial profile
				timePoints      ( 0 3600 31539600 );
				
				// List of axial locations for axial profile
				axialLocations  ( 0.0 0.5 1.0 );
				
				// 2D table for axial profile (z is x axis; time is y axis)
				// The user provides a 2D table (z-vs-time) of normalized fast flux
				// values (normalized to the average).
				// Other options might be available (see axialProfile class).
				data
				(
					( 1.0 1.0 1.0 )
					( 0.5 1.0 1.5 )
					( 0.8 1.1 1.0 )
				);
				
				// Select the axial interpolation method
				axialInterpolationMethod    linear;
				
				// Select the time interpolation method for time steps that fall in between
				// the time points indicated in the timePoints list above
				timeInterpolationMethod     linear;
			}
		
		}
	\endverbatim
	
	As the time point or data lists might be quite long for a realistic history,
	making the solver dictionary less readable, the user can create separate files
	and include them in the fastFluxOptions dictionary using the #include command.
	In the example below the files are stored in the folder `constant/lists`.
	E.g.:
	
	\verbatim
		fastFluxOptions
		{
			timePoints  #include "$FOAM_CASE/constant/lists/timePoints";
			fastFlux    #include "$FOAM_CASE/constant/lists/avgPower";
			timeInterpolationMethod linear;
			
			...
		}
	\endverbatim

SourceFiles
    timeDependentAxialProfile.C

\todo
    - material zones provided or decided based on material.isFuel() or everywhere
    - do we need also thermalFlux classes?
    - do we need the profile to be normalized?
    - do we need a radial profile for the fast/thermal flux?    

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef timeDependentAxialProfile_H
#define timeDependentAxialProfile_H

#include "constantFastFlux.H"
#include "axialProfile.H"
#include "radialProfile.H"
#include "volFields.H"
#include "InterpolateTables.H"
#include "globalOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class timeDependentAxialProfile Declaration
\*---------------------------------------------------------------------------*/

class timeDependentAxialProfile
:
    public constantFastFlux
{
    // Private typedefs
    
        typedef interpolateTableBase::interpolationMethod interpolationMethod;
 
    // Private data
        
        //- Tabulated time-dependent fast flux (n/cm2/s)
        scalarField fluxData_;

        //- Tabulated time values
        scalarField timeData_;
        
        interpolationMethod method_;

        //- Time-dependent rod-average fast flux (n/cm2/s)
        scalarInterpolateTable fluxTable_;
        
        //- Pin axial direction
        vector pinDirection_;
        
        //- Minimum z
        scalar zMin_;
        
        //- Maximum z
        scalar zMax_;
        
        //- Time-dependent axial profile
        autoPtr<axialProfile> axialModel_;

        //- Addressing for fast flux cells
        labelList addr_;
    
    // Private Member Functions

        //- Get addressing for fast flux cells
        void calcAddressing(const dictionary& dict);
        
        //- Calculate the fuel-average cross-sectional area
        void calcMinMaxZ(const dictionary& dict);
        
        //- Disallow default bitwise copy construct
        timeDependentAxialProfile(const timeDependentAxialProfile&);

        //- Disallow default bitwise assignment
        void operator=(const timeDependentAxialProfile&);

public:

    //- Runtime type information
    TypeName("timeDependentAxialProfile");

    // Constructors

        //- Construct from mesh, materials and dict
        timeDependentAxialProfile
        (
            const fvMesh& mesh,
            const materials& materials,
            const dictionary& fastFluxDict
        );


    //- Destructor
    ~timeDependentAxialProfile();


    // Member Functions
        
        //- Update the fastFlux distribution in the supplied list of cells
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
