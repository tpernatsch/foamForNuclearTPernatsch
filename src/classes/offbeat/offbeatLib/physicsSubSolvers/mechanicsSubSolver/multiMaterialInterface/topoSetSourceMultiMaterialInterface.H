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
    Foam::topoSetSourceMultiMaterialInterface

Description
    Interfaces between materials are defined by a topoSetSource,
    e.g. boxToFace, faceToFace.

SourceFiles
    topoSetSourceMultiMaterialInterface.C

\*---------------------------------------------------------------------------*/

#ifndef topoSetSourceMultiMaterialInterface_H
#define topoSetSourceMultiMaterialInterface_H

#include "uniformMultiMaterialInterface.H"
#include "faceSet.H"
#include "topoSetSource.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class topoSetSourceMultiMaterialInterface Declaration
\*---------------------------------------------------------------------------*/

class topoSetSourceMultiMaterialInterface
:
    public uniformMultiMaterialInterface
{
    // Private data
    
        //- Weights to apply at the topoSetSource faces
        scalar interfaceWeight_;
        
        //- List of faces
        faceSet faces_;
        
        //- Source for faceSet
        autoPtr<topoSetSource> source_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        topoSetSourceMultiMaterialInterface(const topoSetSourceMultiMaterialInterface&);

        //- Disallow default bitwise assignment
        void operator=(const topoSetSourceMultiMaterialInterface&);


public:

    //- Runtime type information
    TypeName("topoSetSource");


    // Constructors

        //- Construct from dictionary
        topoSetSourceMultiMaterialInterface
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    ~topoSetSourceMultiMaterialInterface();


    // Member Functions

        //- Correct the face weights
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
