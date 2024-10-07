/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::uniform2DMultiMaterialInterface

Description
    Interfaces between materials are defined by a uniform weight.
    The weight is set to 0 along the z-direction.

SourceFiles
    uniform2DMultiMaterialInterface.C

\*---------------------------------------------------------------------------*/

#ifndef uniform2DMultiMaterialInterface_H
#define uniform2DMultiMaterialInterface_H

#include "uniformMultiMaterialInterface.H"
#include "globalOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class uniform2DMultiMaterialInterface Declaration
\*---------------------------------------------------------------------------*/

class uniform2DMultiMaterialInterface
:
    public uniformMultiMaterialInterface
{
    // Private data

        //- Main axial pin direction
        vector pinDirection_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        uniform2DMultiMaterialInterface(const uniform2DMultiMaterialInterface&);

        //- Disallow default bitwise assignment
        void operator=(const uniform2DMultiMaterialInterface&);


public:

    //- Runtime type information
    TypeName("uniform2D");


    // Constructors

        //- Construct from Dictionary
        uniform2DMultiMaterialInterface
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    ~uniform2DMultiMaterialInterface();


    // Member Functions

        //- Correct the face weights
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
