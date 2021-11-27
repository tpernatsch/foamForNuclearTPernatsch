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

#include "FSPair.H"
#include "TONBModel.H"
#include "myOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TONBModel, 0);
    defineRunTimeSelectionTable
    (
        TONBModel, 
        TONBModels
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TONBModel::TONBModel
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    IOdictionary
    (
        IOobject
        (
            "TONBModel."+typeName,
            pair.mesh().time().timeName(),
            objReg,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    pair_(pair),
    Tsat_(pair.mesh().lookupObject<volScalarField>("T.interface")),
    otherFluidPtr_(nullptr),
    LPtr_(nullptr)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::TONBModel::setPtrs() const
{
    if (otherFluidPtr_ == nullptr)
    {
        HashTable<const fluid*> fluids(pair_.mesh().lookupClass<fluid>());
        otherFluidPtr_ = 
            (fluids[fluids.toc()[0]]->name() == pair_.fluidRef().name()) ?
            fluids[fluids.toc()[1]] : fluids[fluids.toc()[0]];
    }
    if (LPtr_ == nullptr)
    {
        word pairName1(pair_.fluidRef().name()+"."+otherFluidPtr_->name());
        word pairName2(otherFluidPtr_->name()+"."+pair_.fluidRef().name());
        if (pair_.mesh().foundObject<volScalarField>("L."+pairName1))
        {
            LPtr_ = 
                &(
                    pair_.mesh().lookupObject<volScalarField>
                    (
                        "L."+pairName1
                    )
                );
        }
        else
        {
            LPtr_ = 
                &(
                    pair_.mesh().lookupObject<volScalarField>
                    (
                        "L."+pairName2
                    )
                );
        }
    }
}

// ************************************************************************* //
