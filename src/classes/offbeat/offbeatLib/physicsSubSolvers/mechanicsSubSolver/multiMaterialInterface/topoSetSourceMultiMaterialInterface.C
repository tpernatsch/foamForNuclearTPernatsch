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

\*---------------------------------------------------------------------------*/

#include "topoSetSourceMultiMaterialInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(topoSetSourceMultiMaterialInterface, 0);
    addToRunTimeSelectionTable
    (
        multiMaterialInterface, 
        topoSetSourceMultiMaterialInterface, 
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::topoSetSourceMultiMaterialInterface::topoSetSourceMultiMaterialInterface
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    uniformMultiMaterialInterface(mesh, dict),
    interfaceWeight_(readScalar(dict.lookup("interfaceWeights"))),
    faces_
    (
        mesh, 
        "multiMaterialInterface", 
        IOobject::NO_READ, 
        IOobject::NO_WRITE
    ),
    source_()
{
    Istream& is(dict.lookup("source"));
    
    source_ = topoSetSource::New
    (
        word(is),
        mesh,
        is
    );
    
    source_->applyToSet(topoSetSource::ADD, faces_);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::topoSetSourceMultiMaterialInterface::~topoSetSourceMultiMaterialInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- Return a surfaceScalarField of interface weights
void Foam::topoSetSourceMultiMaterialInterface::correct()
{
    uniformMultiMaterialInterface::correct();
    
    //- Default weights without interface
    scalarField& w_i = weights_.ref();
    
    forAllConstIter(faceSet, faces_, iter)
    {
        const label faceI = iter.key();
        
        if (mesh_.isInternalFace(faceI))
        {
            w_i[faceI] = interfaceWeight_;
        }
        else
        {
            NotImplemented("Foam::topoSetSourceMultiMaterialInterface::weights for boundary faces");
        }
    }
}

// ************************************************************************* //
