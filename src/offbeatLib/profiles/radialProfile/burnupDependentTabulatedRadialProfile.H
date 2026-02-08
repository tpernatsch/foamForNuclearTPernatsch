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
    Foam::burnupDependentTabulatedRadialProfile

Description
    Class providing a burnup-dependent radial profile, which is defined by the
    user by interpolation in a 2D lookup table as a function of radial
    location and burnup, i.e. g(Bu,r).
    Both "step" and "linear" interpolation methods are supported (see 
    Foam::InterpolateTable).
    The supplied radial profiles must be normalised to 1, taking note that the
    integral will be in cylindrical coordinates.
    The fuel average burnup is assumed, i.e. no axial effects are included.
    
    The class is provided with a mesh and addressing array, and then deals with
    the mapping of the profile onto an addressed list, which can be accessed
    through the profile function.

Usage
    \verbatim
    {
        type            burnupDependentTabulated;
        burnupPoints    ( 0 70e3 );
        radialLocations ( 0 0.4472135955 0.632455532 0.7745966692 0.894427191 1 );
        data            ( 
            ( 0.8319626974 0.90037028 0.9485969551 1.006237685 1.090645345 1.276336774 )
            ( 0.7179025511 0.8243673183 0.9031851172 1.001405018 1.153001744 1.518179055 )
        );
        radialInterpolationMethod    linear;
        burnupInterpolationMethod    linear;
    }
    \endverbatim

SourceFiles
    burnupDependentTabulatedRadialProfile.C

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    A. Scolaro, E. Brunetto, C. Fiorina - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE 
    LAUSANNE, Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef burnupDependentTabulatedRadialProfile_H
#define burnupDependentTabulatedRadialProfile_H

#include "radialProfile.H"
#include "InterpolateTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class burnupDependentTabulatedRadialProfile Declaration
\*---------------------------------------------------------------------------*/

class burnupDependentTabulatedRadialProfile
:
    public radialProfile
{
    // Private Data
    
        //- Tabulated radial locations, normalised between 0 and 1
        scalarField rValues_;
        
        //- Tabulated burnup values
        scalarField burnupValues_;
        
        //- Tabulated Profile data
        FieldField<Field, scalar> profileData_;
    
        //- Axial location interpolation method
        interpolateTableBase::interpolationMethod rMethod_;
        
        //- Time interpolation method
        interpolateTableBase::interpolationMethod burnupMethod_;
        
        //- Interpolation table
        scalarFieldInterpolateTable table_;
        
    // Private Member Functions

        //- Disallow default bitwise copy construct
        burnupDependentTabulatedRadialProfile(const burnupDependentTabulatedRadialProfile&);

        //- Disallow default bitwise assignment
        void operator=(const burnupDependentTabulatedRadialProfile&);

public:

    //- Runtime type information
    TypeName("burnupDependentTabulated");

    
    // Constructors

        //- Construct from dictionary
        burnupDependentTabulatedRadialProfile
        (
            const fvMesh& mesh,    
            const dictionary& dict,
            const labelList& addr,
            const vector& axialDirection,
            scalar rMax
        );

        
    //- Destructor
    virtual ~burnupDependentTabulatedRadialProfile();

    // Member Functions
    
        //- Update the profile
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
