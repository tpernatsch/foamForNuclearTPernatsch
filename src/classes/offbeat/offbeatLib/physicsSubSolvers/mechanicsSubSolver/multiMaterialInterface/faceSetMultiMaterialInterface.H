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
    Foam::faceSetMultiMaterialInterface

Description
    Interfaces between materials are defined by a faceSet.

SourceFiles
    faceSetMultiMaterialInterface.C

\*---------------------------------------------------------------------------*/

#ifndef faceSetMultiMaterialInterface_H
#define faceSetMultiMaterialInterface_H

#include "uniformMultiMaterialInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class faceSetMultiMaterialInterface Declaration
\*---------------------------------------------------------------------------*/

class faceSetMultiMaterialInterface
:
    public uniformMultiMaterialInterface
{
    // Private data
    
        //- Weights to apply at the faceSet faces
        scalar interfaceWeight_;
        
        //- Reference to faceSet
        faceSet faces_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        faceSetMultiMaterialInterface(const faceSetMultiMaterialInterface&);

        //- Disallow default bitwise assignment
        void operator=(const faceSetMultiMaterialInterface&);


public:

    //- Runtime type information
    TypeName("faceSet");


    // Constructors

        //- Construct from dictionary
        faceSetMultiMaterialInterface
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    ~faceSetMultiMaterialInterface();


    // Member Functions

        //- Correct the face weights
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
