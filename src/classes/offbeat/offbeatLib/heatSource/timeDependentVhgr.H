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
    Volumetric heat source is provided

    \htmlinclude heatSource_timeDependentVhgr.html


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

#include "heatSource.H"
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
    public heatSource
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
            const dictionary& heatSourceOptDict
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
