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
    Foam::constantTabulatedAzimuthalProfile

Description
    Azimuthal profile defined using an constant interpolation table as a
    function of the azimuthal position.
    Both "step" and "linear" interpolation methods are supported (see 
    Foam::InterpolateTable).
    The supplied profile must be normalised to 1.

Usage
    \verbatim
    {
        type            constantTabulated;
        azimuthalLocations ( 0 0.4472135955 0.632455532 0.7745966692 0.894427191 1 );
        data            ( 0.8319626974 0.90037028 0.9485969551 1.006237685 1.090645345 1.276336774 );
        azimuthalInterpolationMethod    linear;
    }
    \endverbatim

SourceFiles
    constantTabulatedAzimuthalProfile.C

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    A. Scolaro, E. Brunetto, C. Fiorina - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE 
    LAUSANNE, Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef constantTabulatedAzimuthalProfile_H
#define constantTabulatedAzimuthalProfile_H

#include "azimuthalProfile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class constantTabulatedAzimuthalProfile Declaration
\*---------------------------------------------------------------------------*/

class constantTabulatedAzimuthalProfile
:
    public azimuthalProfile
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        constantTabulatedAzimuthalProfile(const constantTabulatedAzimuthalProfile&);

        //- Disallow default bitwise assignment
        void operator=(const constantTabulatedAzimuthalProfile&);

public:

    //- Runtime type information
    TypeName("constantTabulated");

    
    // Constructors

        //- Construct from dictionary
        constantTabulatedAzimuthalProfile
        (
            const fvMesh& mesh,    
            const dictionary& dict,
            const labelList& addr,
            const vector& axialDirection,
            const vector& radialDirection
        );

        
    //- Destructor
    virtual ~constantTabulatedAzimuthalProfile();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
