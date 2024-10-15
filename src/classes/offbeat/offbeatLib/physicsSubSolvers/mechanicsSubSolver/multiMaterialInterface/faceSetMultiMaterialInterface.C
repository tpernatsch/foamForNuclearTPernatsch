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

#include "faceSetMultiMaterialInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceSetMultiMaterialInterface, 0);
    addToRunTimeSelectionTable
    (
        multiMaterialInterface, 
        faceSetMultiMaterialInterface, 
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceSetMultiMaterialInterface::faceSetMultiMaterialInterface
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    uniformMultiMaterialInterface(mesh, dict),
    interfaceWeight_(readScalar(dict.lookup("interfaceWeights"))),
#ifdef OPENFOAMFOUNDATION
    faces_(mesh, dict.lookup("faceSetName"))
#elif OPENFOAMESI
    faces_(mesh, dict.get<word>("faceSetName"))
#endif    
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceSetMultiMaterialInterface::~faceSetMultiMaterialInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- Return a surfaceScalarField of interface weights
void Foam::faceSetMultiMaterialInterface::correct()
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
            NotImplemented("Foam::faceSetMultiMaterialInterface::weights for boundary faces");
        }
    }
}

// ************************************************************************* //
