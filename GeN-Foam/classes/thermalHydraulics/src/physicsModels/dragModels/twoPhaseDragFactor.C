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

#include "twoPhaseDragFactor.H"
#include "myOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseDragFactor, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseDragFactor::twoPhaseDragFactor
(
    FSPair& F1SPair,
    FSPair& F2SPair,
    const dictionary& dict
)
:
    IOdictionary
    (
        IOobject
        (
            typeName,
            F1SPair.mesh().time().timeName(),
            F1SPair.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(F1SPair.mesh()),
    F1SPair_(F1SPair),
    F2SPair_(F2SPair),
    KdTotU_
    (
        IOobject
        (
            "KdTotU",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor
        (
            "", 
            dimDensity*dimVelocity/dimTime, 
            tensor
            (
                1, 0, 0,
                0, 1, 0,
                0, 0, 1
            )
        ),
        zeroGradientFvPatchTensorField::typeName
    ),
    twoPhaseDragMultiplierPtr_
    (
        twoPhaseDragMultiplierModel::New
        (
            mesh_,
            dict.subDict
            (
                "physicsModels"
            ).subDict
            (
                "twoPhaseDragMultiplierModel"
            ),
            mesh_
        )
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoPhaseDragFactor::correct()
{
    bool notFirstTimeStep
    (
        (mesh_.time().timeIndex() != mesh_.time().startTimeIndex()+1)
    );
    if (notFirstTimeStep)
        myOps::storePrevIterIfRelax(KdTotU_);
    twoPhaseDragMultiplierPtr_->correctField(KdTotU_);
    if (notFirstTimeStep)
        KdTotU_.relax();
    dimensionedScalar minMagU("", dimVelocity, 1e-9);
    F1SPair_.Kd() = 
        F1SPair_.f()*KdTotU_/max(F1SPair_.fluidRef().magU(), minMagU);
    F2SPair_.Kd() = 
        F2SPair_.f()*KdTotU_/max(F2SPair_.fluidRef().magU(), minMagU);
}

// ************************************************************************* //
