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
    Foam::timeDependentLhgr

Description
	Similar to the `constantLhgr` heatSource class, the `timeDependentLhgr` class in
	OFFBEAT allows you to set a average **linear heat generation rate (lhgr)**.
	The main difference is that the lhgr provided by the user can vary over time.
	It is also possible to apply a radial and an axial profile to the lhgr. Note
	that even if the lhgr is constant over time, the axial and radial profile might
	change over time, depending on the profile type chosen by the user.
	
	#### Formulation
	
	The final volumetric heat source is derived as a combination of lhgr, radial and
	axial profiles, as follows:

	$$
	\begin{aligned}
	Q(t,r,z) = Q_{avg}(t) \cdot f(r,t) \cdot g(z,t)
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
	Q_{avg}(t) = \frac{lhgr(t) \cdot h_{ref}}{\frac{V_{model}}{\alpha}}
	\end{aligned}
	$$
	
	where:
	
	- \f$h_{ref}\f$ is the height of the rod in the reference, undeformed geometrY
	- \f$V_{model}\f$ is the volume of the fuel as calculated by OpenFOAM in the
	model
	- \f$\alpha\f$ is the fraction of the real geometry modeled in the simulations
	(e.g. 0.25 if you are using symmetry BCs and you are simulating only a quarter
	of the rod).
	
________________________________________________________________________________

\par Options
	The `timeDependentLhgr` heat source model requires the user to specify a few
	additional parameters in the `heatSourceOptions` sub-dictionary located 
	within the `solverDict`.
    
    Parameters in `heatSourceOptions`:
	
	- <b>`timePoints`</b> - A list of time values at which the lhgr is provided. The time
	unit depend on the userTime selected by the user (seconds by default).
	- <b>`lhgr`</b> - A list of lhgr values in W/m, one value per time-point.
	- <b>`timeInterpolationMethod`</b> - Select the time interpolation method for time
	steps that fall in between the time points indicated in the timePoints list.
	- <b>`materials`</b> - A list of material (or cellZones names) where the heat source
	model applies.
	- <b>`axialProfile`</b> - A sub-dictionary that specifies the type of axial profile.
	- <b>`radialProfile`</b> - A sub-dictionary that specifies the type of radial profile.
	- <b>`useReferenceDimensions`</b> – A bool to activate reading reference dimensions for 
	axial profile (`zMin` and `zMax`) from the dictionary (false by default).  
	- <b>`zMin`</b> – The reference axial coordinate of the bottom of the fuel region. 
	Only used if `useReferenceDimensions` is set to true.  
	- <b>`zMax`</b> – The reference axial coordinate of the top of the fuel region. 
	Only used if `useReferenceDimensions` is set to true.
	- <b>`volumeFraction`</b> – The fraction of the reference rod present in the local 
	mesh. Only used if `useReferenceDimensions` is set to true.
	
	The `timeDependentLhgr` heat source model is designed to calculate the
	volumetric power density for a cylindrical rod with an initial reference height
	(i.e., in undeformed geometry) along a given axial direction. For simulations
	where there is an axis of symmetry, the class can take into consideration the
	angular fraction of the real rod that is being modeled.
	
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
	Here is a code snippet of the `solverDict` to be used for the case of a lhgr
	that ramps to 10 kW/m in 1hr time, remains constant for a year, and then ramps
	down to 0 in 1 hr.
	This particular example considers  flat (i.e. uniform) axial and radial
	profiles, and a quarter of a fuel disc model.
	
	\verbatim
		heatSource timeDependentLhgr;

		globalOptions
		{
			angularFraction 0.25;
			pinDirection  (0 0 1);
		}

		heatSourceOptions
		{
			//           t0  1hr    1yr       1yr+1hr
			timePoints  ( 0   3600   31536000  31539600);
			lhgr        ( 0   100e2  100e2     0       );

			timeInterpolationMethod linear; //step;

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
				// Other options might be available (see axialProfile class).
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

				// Select the time interpolation method for time steps that fall in between
				// the time points indicated in the timePoints list above
				timeInterpolationMethod     linear;
			}
		}
	\endverbatim
	
	As the time point or data lists might be quite long for a realistic history,
	making the solver dictionary less readable, the user can create separate files
	and include them in the heatSourceOptions dictionary using the #include command.
	E.g.:
	
	\verbatim
		heatSourceOptions
		{
			timePoints  #include "$FOAM_CASE/constant/lists/timePoints";
			lhgr        #include "$FOAM_CASE/constant/lists/avgPower";
			timeInterpolationMethod linear;

			...
		}
	\endverbatim

SourceFiles
    timeDependentLhgr.C

\todo
    - The conversion from lhgr to volumetric power might be wrong when fuel with 
    different dimension (e.g. hollow and full pellet in the same rod) are used

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef timeDependentLhgr_H
#define timeDependentLhgr_H

#include "constantLhgr.H"
#include "axialProfile.H"
#include "radialProfile.H"
#include "volFields.H"
#include "InterpolateTables.H"
#include "globalOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class timeDependentLhgr Declaration
\*---------------------------------------------------------------------------*/

class timeDependentLhgr
:
    public constantLhgr
{
    // Private typedefs
    
        typedef interpolateTableBase::interpolationMethod interpolationMethod;
 
    // Private data
        
        //- Tabulated time-dependent rod-average LHGR (W/m)
        scalarField lhgrData_;

        //- Tabulated time values
        scalarField timeData_;
        
        interpolationMethod method_;

        //- Time-dependent table for rod-average LHGR (W/m)
        scalarInterpolateTable lhgrTable_;

        //- Maximum allowed Lhgr increase for adjustable time step
        scalar maxAbsoluteLhgrIncrease_;
        
    // Private Member Functions

        //- Lookup the LHGR from the supplied physical time in seconds
        inline scalar lookupLHGR(scalar t) const
        {
            return lhgrTable_(mesh_.time().timeToUserTime(t));
        }
        
        //- Disallow default bitwise copy construct
        timeDependentLhgr(const timeDependentLhgr&);

        //- Disallow default bitwise assignment
        void operator=(const timeDependentLhgr&);

public:

    //- Runtime type information
    TypeName("timeDependentLhgr");

    // Constructors

        //- Construct from mesh, materials and dict
        timeDependentLhgr
        (
            const fvMesh& mesh,
            const materials& materials,
            const dictionary& heatSourceDict
        );


    //- Destructor
    ~timeDependentLhgr();


    // Member Functions
    
        //- update the power distribution in the supplied list of cells
        virtual void correct();
        
        //- Return the next average power density (assuming the deltaT stays the
        // same)
        virtual scalar predictNextPowerDensity() const;    

        //- Provide next deltaT due to absolute change in heat source (not relative)
        virtual scalar nextDeltaTAbsoluteChange(scalar deltaT) const;   

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
