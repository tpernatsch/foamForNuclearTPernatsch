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
    Foam::timeDependentVhgr

Description
	The `timeDependentVhgr` class in OFFBEAT allows you to set a average
	**volumetric heat generation rate (vhgr)** in unit of <code>W/m<sup>3</sup></code>. 
    Similar	to `timeDependentVhgr`, this class allow users to provide time-dependent vhgr.
	Since this class is mainly developed for TRISO, it currently doesn't support the
	addition of a radial and an axial profile of vhgr (as done in the lhgr heat
	source classes).	
	
________________________________________________________________________________

\par Options	
	The `timeDependentVhgr` heat source model requires the user to specify a few
	additional parameters in the `heatSourceOptions` sub-dictionary located 
	within the `solverDict`.
    
    Parameters in `heatSourceOptions`:
	
	- <b>`timePoints`</b> - A list of time values at which the vhgr is provided. 
    The time unit depend on the userTime selected by the user (seconds by default).
	- <b>`vhgr`</b> - A list of vhgr values in <code>W/m<sup>3</sup></code>, one 
    value per time-point.
	- <b>`timeInterpolationMethod`</b> - Select the time interpolation method for time
	steps that fall in between the time points indicated in the timePoints list.
	- <b>`materials`</b> - A list of material (or cellZones names) where the heat source
	model applies.
	
________________________________________________________________________________
	
Usage
	Here is a code snippet of the `solverDict` to be used for the case of a vhgr
	that ramps to 1e9 <code>W/m<sup>3</sup></code> in 1hr time, remains constant for 
    a year, and then ramps down to 0 in 1 hr.
	
	\verbatim
        heatSource timeDependentVhgr;
	
        heatSourceOptions
        {
            //           t0  1hr    1yr       1yr+1hr
            timePoints  ( 0   3600 31536000  31539600);
            vhgr        ( 0   1e9  1e9       0       );

            timeInterpolationMethod linear; //step;

            materials ( fuel );
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
    timeDependentVhgr.C

\mainauthor
    F. Xiang - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour) & XJTU (Xi'an JiaoTong
    University, China, Nuclear Thermal-hydraulic Laboratory)
    A. Scolaro - EPFL

\contribution
    E. Brunetto - EPFL

\date
    Feburary 2023

\*---------------------------------------------------------------------------*/

#ifndef timeDependentVhgr_H
#define timeDependentVhgr_H

#include "constantHeatSource.H"
#include "axialProfile.H"
#include "radialProfile.H"
#include "volFields.H"
#include "InterpolateTables.H"
#include "globalOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class timeDependentVhgr Declaration
\*---------------------------------------------------------------------------*/

class timeDependentVhgr
:
    public constantHeatSource
{
    // Private typedefs

        typedef interpolateTableBase::interpolationMethod interpolationMethod;

    // Private data

        //- Tabulated time-dependent  vhgr (W/m³)
        scalarField vhgrData_;

        //- Tabulated time values
        scalarField timeData_;

        interpolationMethod method_;

        //- Time-dependent table for vhgr (W/m³)
        scalarInterpolateTable vhgrTable_;

        //- Addressing for power-producing cells
        labelList addr_;

    // Private Member Functions

        //- Get addressing for power-producing cells
        void calcAddressing(const dictionary& dict);

        //- Lookup the vhgr from the supplied physical time in seconds
        inline scalar lookupVHGR(scalar t) const
        {
            return vhgrTable_(mesh_.time().timeToUserTime(t));
        }

        //- Disallow default bitwise copy construct
        timeDependentVhgr(const timeDependentVhgr&);

        //- Disallow default bitwise assignment
        void operator=(const timeDependentVhgr&);

public:

    //- Runtime type information
    TypeName("timeDependentVhgr");

    // Constructors

        //- Construct from mesh, materials and dict
        timeDependentVhgr
        (
            const fvMesh& mesh,
            const materials& materials,
            const dictionary& heatSourceDict
        );


    //- Destructor
    ~timeDependentVhgr();


    // Member Functions

        //- update the power distribution in the supplied list of cells
        virtual void correct();

        //- Return the next average power density (assuming the deltaT stays the
        // same)
        virtual scalar predictNextPowerDensity() const;

        //- Return the average power density
        virtual scalar averagePowerDensity() const;

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
