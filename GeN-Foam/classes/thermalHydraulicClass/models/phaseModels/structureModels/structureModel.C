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

#include "structureModel.H"
#include "fvmDdt.H"
#include "fvmSup.H"
#include "fluid.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(structureModel, 0);
    defineRunTimeSelectionTable
    (
        structureModel,
        structureModels
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::structureModel::structureModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    phaseBase
    ( 
        dict,
        mesh,
        "structure"
    ),
    regions_(0),
    tortuosity_
    (
        IOobject
        (
            "tortuosity",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedTensor("Dh", dimless, tensor::one),
        zeroGradientFvPatchScalarField::typeName
    ),
    Dh_
    (
        IOobject
        (
            "Dh",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("Dh", dimLength, SMALL),
        zeroGradientFvPatchScalarField::typeName
    ),
    lDh_
    (
        IOobject
        (
            "lDh",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector("lDh", dimLength, vector(SMALL, SMALL, SMALL)),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tact_
    (
        IOobject
        (
            "T.activeStructure",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ, //- Constructed from powerModel initial Ts
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("T.structure", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    iAact_
    (
        IOobject
        (
            "iA.activeStructure",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("iA.structure", dimArea/dimVol, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tpas_
    (
        IOobject
        (
            "T.passiveStructure",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("T.passiveStructure", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    iApas_
    (
        IOobject
        (
            "iA.passiveStructure",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("iA.passiveStructure", dimArea/dimVol, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    rhopas_
    (
        IOobject
        (
            "rho.passiveStructure",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("rho.passiveStructure", dimMass/dimVol, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Cppas_
    (
        IOobject
        (
            "Cp.passiveStructure",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar
        (
            "Cp.passiveStructure", 
            dimEnergy/dimMass/dimTemperature, 
            0
        ),
        zeroGradientFvPatchScalarField::typeName
    ),
    Rg2l_
    (
        IOobject
        (
            "R.global2Local",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedTensor
        (
            "", 
            dimless, 
            tensor
            (
                1,0,0,
                0,1,0,
                0,0,1
            )
        ),
        zeroGradientFvPatchScalarField::typeName
    ),
    Rl2g_
    (
        IOobject
        (
            "R.local2Global",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedTensor
        (
            "", 
            dimless, 
            tensor
            (
                1,0,0,
                0,1,0,
                0,0,1
            )
        ),
        zeroGradientFvPatchScalarField::typeName
    ),
    nuclearFuelPowerDensity_
    (
        IOobject
        (
            "nuclearFuelPowerDensity",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "nuclearFuelPowerDensity", 
            dimPower/dimVol, 
            0
        ),
        zeroGradientFvPatchScalarField::typeName
    ),
    TavFuel_
    (
        IOobject
        (
            "TavFuel",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "TavFuel", 
            dimTemperature, 
            0
        ),
        zeroGradientFvPatchScalarField::typeName
    ),
    TavClad_
    (
        IOobject
        (
            "TavClad",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "TavClad", 
            dimTemperature, 
            0
        ),
        zeroGradientFvPatchScalarField::typeName
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::structureModel::correct
(
    const volScalarField& HT,
    const volScalarField& H
)
{
    //- Correct powerModels & correct active structure surface temperature
    Tact_ *= 0.0;
    forAllIter
    (
        powerModelTable,
        powerModels_,
        iter
    )
    {
        iter()->correct(HT, H);
        Tact_ += iter()->T();
    }
    Tact_.correctBoundaryConditions();

    //- Correct inert subStructure
    volScalarField alpha
    (
        Foam::max(*this, dimensionedScalar("", dimless, 1e-69))
    );
    fvScalarMatrix pasEqn
    (
            fvm::ddt(alpha*rhopas_*Cppas_, Tpas_)
        ==
            iApas_*HT
        -   fvm::Sp(iApas_*H, Tpas_)
    );
    pasEqn.solve();
}

Foam::tmp<Foam::volScalarField> Foam::structureModel::explicitHeatSource
(
    const fluid& fluid,
    const volScalarField& H
) const
{
    const volScalarField& T(fluid.thermo().T());
    tmp<volScalarField> tQ
    (
        new volScalarField
        (
                iAact_*H*(Tact_-T)
            +   iApas_*H*(Tpas_-T)
        )
    );
    volScalarField& Q = tQ.ref();
    Q.correctBoundaryConditions();
    return tQ;
}


Foam::tmp<Foam::fvScalarMatrix> 
Foam::structureModel::linearizedSemiImplicitHeatSource
(
    fluid& fluid,
    const volScalarField& H
) const
{
    volScalarField& he(fluid.thermo().he());
    const volScalarField& T(fluid.thermo().T());
    tmp<volScalarField> Cp(fluid.thermo().Cp());
    tmp<fvScalarMatrix> tQ
    (
        new fvScalarMatrix
        (
                //- Source/sink due to active subStructure
                iAact_*H*(Tact_-T+he/Cp())
            -   fvm::Sp(iAact_*H/Cp(), he)
                //- Source/sink due to passive subStructure
            +   iApas_*H*(Tpas_-T+he/Cp())
            -   fvm::Sp(iApas_*H/Cp(), he)
        )
    );
    return tQ;
}

// ************************************************************************* //
