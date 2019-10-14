/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "constantHeatTransfer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable
    (
        heatTransferModel, 
        constant, 
        FFHeatTransferModels
    );
    addToRunTimeSelectionTable
    (
        heatTransferModel, 
        constant, 
        FSHeatTransferModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::constant::constant
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FFPair& FFPair,
    const word nameBulk,
    const word nameInterface
)
:
    heatTransferModel
    (
        objReg,
        dict,
        FFPair,
        nameBulk,
        nameInterface
    ),
    htc_
    (
        IOobject
        (
            IOobject::groupName(typeName, this->pairName_),
            FFPair.mesh().time().timeName(),
            objReg
        ),
        FFPair.mesh(),
        dimensionedScalar("value", this->dimHtc, dict)
    )
{}

Foam::heatTransferModels::constant::constant
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FSPair& FSPair,
    const word region
)
:
    heatTransferModel
    (
        objReg,
        dict,
        FSPair,
        region
    ),
    htc_
    (
        IOobject
        (
            IOobject::groupName(typeName, this->pairName_),
            FSPair.mesh().time().timeName(),
            objReg
        ),
        FSPair.mesh(),
        dimensionedScalar("value", this->dimHtc, dict),
        zeroGradientFvPatchScalarField::typeName
    )
{
    //- Correct htc so that it is 0 outside the bounds of the given structure
    //  region
    htc_.primitiveFieldRef() *= FSPair.structure().cellFields()[region];
    htc_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::heatTransferModels::constant::htc() const
{
    tmp<volScalarField> thtc
    (
        new volScalarField(htc_)
    );

    return thtc;
}


// ************************************************************************* //
