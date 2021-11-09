/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "FFPair.H"
#include "saturationModel.H"
#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(saturationModel, 0);
    defineRunTimeSelectionTable(saturationModel, saturationModels);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModel::saturationModel
(
    const phaseChangeModel& pcm,
    const dictionary& dict, 
    const objectRegistry& objReg
)
:
    IOdictionary
    (
        IOobject
        (
            "saturationModel",
            pcm.mesh().time().constant(),
            pcm.mesh()
        ),
        dict
    ),
    mesh_(pcm.mesh()),
    liquid_(pcm.liquid()),
    vapour_(pcm.vapour())
{
    map_["TSat"] = &Foam::saturationModel::valueTSat;
    map_["pSat"] = &Foam::saturationModel::valuePSat;
    map_["lnPSat"] = &Foam::saturationModel::valueLnPSat;
    map_["pSatPrime"] = &Foam::saturationModel::valuePSatPrime;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::saturationModel::correctField
(
    volScalarField& field, 
    const word& func
)
{
    funcPtr_ = map_[func];
    forAll(mesh_.cells(), i)
    {
        field[i] = (this->*funcPtr_)(i);
    }
    field.correctBoundaryConditions();
}

// ************************************************************************* //
