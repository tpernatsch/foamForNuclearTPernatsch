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
#include "fluidDiameterModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidDiameterModel, 0);
    defineRunTimeSelectionTable
    (
        fluidDiameterModel, 
        fluidDiameterModels
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidDiameterModel::fluidDiameterModel
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
            typeName,
            pair.mesh().time().timeName(),
            objReg,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(pair.mesh()),
    pair_(pair),
    fluid_(pair.fluidRef()),
    dispersion_(pair.fluidRef().dispersion()),
    Dhs_(pair.structureRef().Dh()),
    Dh_(pair.fluidRef().Dh())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidDiameterModel::~fluidDiameterModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluidDiameterModel::correctField(volScalarField& Dh) const
{
    forAll(mesh_.cells(), i)
    {
        const scalar& di(dispersion_[i]);
        if (di == 0)
            Dh[i] = Dhs_[i];
        else
            Dh[i] = di*this->value(i) + (1.0-di)*Dhs_[i];
    }
    Dh.correctBoundaryConditions();
}


// ************************************************************************* //
