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

#include "uniformMultiMaterialInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(uniformMultiMaterialInterface, 0);
    addToRunTimeSelectionTable
    (
        multiMaterialInterface, 
        uniformMultiMaterialInterface, 
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformMultiMaterialInterface::uniformMultiMaterialInterface
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    multiMaterialInterface(mesh, dict),
    defaultWeight_(readScalar(dict.lookup("defaultWeights")))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::uniformMultiMaterialInterface::~uniformMultiMaterialInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- Return a surfaceScalarField of interface weights
void Foam::uniformMultiMaterialInterface::correct()
{
    //- Default internalField weights without interface
    scalarField& w_i = weights_.ref();
    w_i = defaultWeight_;
    
    //- Apply weights to coupled boundaries
    forAll(weights_.boundaryField(), patchI)
    {
        fvsPatchField<scalar>& pf = weights_.boundaryFieldRef()[patchI];
        
        if(pf.coupled())
        {
            pf = defaultWeight_;
        }
        else
        {
            pf = 0.0;
        }
    }
}

// ************************************************************************* //
