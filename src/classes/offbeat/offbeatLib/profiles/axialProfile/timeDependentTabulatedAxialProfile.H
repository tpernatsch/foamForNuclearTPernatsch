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
    Foam::timeDependentTabulatedAxialProfile

Description
    Class providing a time-dependent axial profile, which is defined by the
    user by interpolation in a 2D lookup table as a function of axial
    location and time, i.e. f(z,t).
    Both "step" and "linear" interpolation methods are supported (see 
    Foam::InterpolateTable).
    The supplied profile must be normalised to 1.
    
    The class is provided with a mesh and addressing array, and then deals with
    the mapping of the profile onto an addressed list, which can be accessed
    through the profile function.

Usage
    \verbatim
    {
        type            timeDependentTabulated;
        timePoints      ( 0 100 );
        axialLocations  ( 0.0 0.5 1.0 );
        data          ( ( 0.5 1.0 1.5 )
                        ( 0.8 1.1 1.0 ) );
        axialInterpolationMethod    linear;
        timeInterpolationMethod     linear;
        
    }
    \endverbatim

SourceFiles
    timeDependentTabulatedAxialProfile.C

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    A. Scolaro, E. Brunetto, C. Fiorina - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE 
    LAUSANNE, Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef timeDependentTabulatedAxialProfile_H
#define timeDependentTabulatedAxialProfile_H

#include "axialProfile.H"
#include "InterpolateTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class timeDependentTabulatedAxialProfile Declaration
\*---------------------------------------------------------------------------*/

class timeDependentTabulatedAxialProfile
:
    public axialProfile
{
    // Private Data
    
        //- Tabulated axial locations, normalised between 0 and 1
        scalarField zValues_;
        
        //- Tabulated time values
        scalarField timeValues_;
        
        //- Tabulated Profile data
        FieldField<Field, scalar> profileData_;
    
        //- Axial location interpolation method
        interpolateTableBase::interpolationMethod zMethod_;
        
        //- Time interpolation method
        interpolateTableBase::interpolationMethod tMethod_;
        
        //- Interpolation table
        scalarFieldInterpolateTable table_;
        
        
    // Private Member Functions

        //- Disallow default bitwise copy construct
        timeDependentTabulatedAxialProfile(const timeDependentTabulatedAxialProfile&);

        //- Disallow default bitwise assignment
        void operator=(const timeDependentTabulatedAxialProfile&);

public:

    //- Runtime type information
    TypeName("timeDependentTabulated");

    
    // Constructors

        //- Construct from dictionary
        timeDependentTabulatedAxialProfile
        (
            const fvMesh& mesh,    
            const dictionary& dict,
            const labelList& addr,
            const vector& axialDirection,
            scalar zMin,
            scalar zMax
        );

        
    //- Destructor
    virtual ~timeDependentTabulatedAxialProfile();

    // Member Functions
    
        //- Update the profile
        virtual void correct();
        
        //- For time-dependent heat sources, return the next time marker,
        //- i.e. the next time point at which the rate of change
        //- will be discontinuous or non-smooth.
        virtual scalar nextTimeMarker() const
        {
            return mesh_.time().userTimeToTime(
                table_.getNextPoint(
                    mesh_.time().timeToUserTime(mesh_.time().value())
                )
            );
        }

        //- Return the final time marker after which the field distribution
        //- is no longer defined.
        virtual scalar lastTimeMarker() const
        {
            return mesh_.time().userTimeToTime(table_.getLastPoint());
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
