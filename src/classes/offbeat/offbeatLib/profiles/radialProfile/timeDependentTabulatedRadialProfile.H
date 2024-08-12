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
    Foam::timeDependentTabulatedRadialProfile

Description
    Class providing a time-dependent radial profile, which is defined by the
    user by interpolation in a 2D lookup table as a function of radial
    location and time, i.e. g(r,t).
    Both "step" and "linear" interpolation methods are supported (see 
    Foam::InterpolateTable).
    The supplied profiles must be normalised to 1, taking note that the
    integral is in cylindrical coordinates.
    
    The class is provided with a mesh and addressing array, and then deals with
    the mapping of the profile onto an addressed list, which can be accessed
    through the profile function.

Usage
    \verbatim
    {
        type            timeDependentTabulated;
        timePoints      ( 0 100 );
        radialLocations ( 0 0.4472135955 0.632455532 0.7745966692 0.894427191 1 );
        data            ( 
            ( 0.8319626974 0.90037028 0.9485969551 1.006237685 1.090645345 1.276336774 )
            ( 0.7179025511 0.8243673183 0.9031851172 1.001405018 1.153001744 1.518179055 )
        );
        radialInterpolationMethod    linear;
        timeInterpolationMethod     linear;
    }
    \endverbatim

SourceFiles
    timeDependentTabulatedRadialProfile.C

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    A. Scolaro, E. Brunetto, C. Fiorina - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE 
    LAUSANNE, Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef timeDependentTabulatedRadialProfile_H
#define timeDependentTabulatedRadialProfile_H

#include "radialProfile.H"
#include "InterpolateTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class timeDependentTabulatedRadialProfile Declaration
\*---------------------------------------------------------------------------*/

class timeDependentTabulatedRadialProfile
:
    public radialProfile
{
    // Private Data
    
        //- Tabulated radial locations, normalised between 0 and 1
        scalarField rValues_;
        
        //- Tabulated time values
        scalarField timeValues_;
        
        //- Tabulated Profile data
        FieldField<Field, scalar> profileData_;
    
        //- Axial location interpolation method
        interpolateTableBase::interpolationMethod rMethod_;
        
        //- Time interpolation method
        interpolateTableBase::interpolationMethod tMethod_;
        
        //- Interpolation table
        scalarFieldInterpolateTable table_;
        
        
    // Private Member Functions

        //- Disallow default bitwise copy construct
        timeDependentTabulatedRadialProfile(const timeDependentTabulatedRadialProfile&);

        //- Disallow default bitwise assignment
        void operator=(const timeDependentTabulatedRadialProfile&);

public:

    //- Runtime type information
    TypeName("timeDependentTabulated");

    
    // Constructors

        //- Construct from dictionary
        timeDependentTabulatedRadialProfile
        (
            const fvMesh& mesh,    
            const dictionary& dict,
            const labelList& addr,
            const vector& axialDirection,
            scalar rMax
        );

        
    //- Destructor
    virtual ~timeDependentTabulatedRadialProfile();

    // Member Functions
    
        //- Update the profile
        virtual void correct();

        //- For time-dependent heat sources, return the next time marker,
        //- i.e. the next time point at which the rate of power change
        //- will be discontinuous or non-smooth.
        virtual scalar nextTimeMarker() const
        {
            return mesh_.time().userTimeToTime(
                table_.getNextPoint(
                    mesh_.time().timeToUserTime(mesh_.time().value())
                )
            );
        }

        //- Return the final time marker after which the power distribution
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
