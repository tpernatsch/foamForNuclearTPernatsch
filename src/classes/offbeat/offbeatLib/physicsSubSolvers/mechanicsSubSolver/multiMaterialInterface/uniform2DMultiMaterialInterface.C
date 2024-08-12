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

#include "uniform2DMultiMaterialInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(uniform2DMultiMaterialInterface, 0);
    addToRunTimeSelectionTable
    (
        multiMaterialInterface, 
        uniform2DMultiMaterialInterface, 
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniform2DMultiMaterialInterface::uniform2DMultiMaterialInterface
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    uniformMultiMaterialInterface(mesh, dict),
    pinDirection_(mesh_.lookupObject<globalOptions>("globalOptions").pinDirection())
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::uniform2DMultiMaterialInterface::~uniform2DMultiMaterialInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- Return a surfaceScalarField of interface weights
void Foam::uniform2DMultiMaterialInterface::correct()
{
    Foam::uniformMultiMaterialInterface::correct();

    const fvMesh& referenceMesh = 
    ( 
        mesh_.foundObject<fvMesh>("referenceMesh")
    ) ?
    mesh_.lookupObject<fvMesh>("referenceMesh")
    : 
    mesh_;

    scalarField& w_i = weights_.ref();

    surfaceVectorField nf = referenceMesh.Sf()/referenceMesh.magSf();
    
    forAll(w_i, faceI)
    {
        if ( mag(((nf[faceI] & pinDirection_) - 1.0)) < SMALL )
        {
            w_i[faceI] = 0 ;
        }
    }
}

// ************************************************************************* //
