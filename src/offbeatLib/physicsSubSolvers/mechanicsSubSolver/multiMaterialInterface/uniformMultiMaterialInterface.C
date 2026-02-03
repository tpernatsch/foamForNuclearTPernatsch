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
    defaultWeightEps_(readScalar(dict.lookup("defaultWeights"))),
    defaultWeightGrad_(defaultWeightEps_)
{
    if (dict.found("defaultWeightsGrad"))
    {
        defaultWeightGrad_ = readScalar(dict.lookup("defaultWeightsGrad"));
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::uniformMultiMaterialInterface::~uniformMultiMaterialInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- Return a surfaceScalarField of interface weights
void Foam::uniformMultiMaterialInterface::correct()
{
    //- Default internalField weights without interface
    scalarField& wEps_i = weightsEps_.ref();
    scalarField& wGrad_i = weightsGrad_.ref();
    wEps_i = defaultWeightEps_;
    wGrad_i = defaultWeightGrad_;
    
    //- Apply weights to coupled boundaries
    forAll(weightsEps_.boundaryField(), patchI)
    {
        fvsPatchField<scalar>& pfEps = weightsEps_.boundaryFieldRef()[patchI];
        fvsPatchField<scalar>& pfGrad = weightsGrad_.boundaryFieldRef()[patchI];
        
        if(pfEps.coupled())
        {
            pfEps = defaultWeightEps_;
            pfGrad = defaultWeightGrad_;
        }
        else
        {
            pfEps = defaultWeightEps_;
            pfGrad = defaultWeightGrad_;
        }
    }
}

// ************************************************************************* //
