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
    cells_(0),
    tortuosity_
    (
        IOobject
        (
            "tortuosity",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedTensor("", dimless, tensor::one),
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
        dimensionedScalar("", dimLength, SMALL),
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
        dimensionedVector("", dimLength, vector(SMALL, SMALL, SMALL)),
        zeroGradientFvPatchScalarField::typeName
    ),
    heatFlux_
    (
        IOobject
        (
            "heatFlux.structure",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimPower/dimArea, 0.0),
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
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    iAact_
    (
        IOobject
        (
            "volumetricArea.activeStructure",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("", dimArea/dimVol, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    alphapas_
    (
        IOobject
        (
            "alpha.passiveStructure",
            mesh.time().timeName(),
            mesh
        ),
        *this
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
        dimensionedScalar("T", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    iApas_
    (
        IOobject
        (
            "volumetricArea.passiveStructure",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("", dimArea/dimVol, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    alphaRhoCppas_
    (
        IOobject
        (
            "alphaRhoCp.passiveStructure",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar
        ("", dimEnergy/dimVol/dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Twall_
    (
        IOobject
        (
            "T.wall",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("T", dimTemperature, 0),
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
    )
{
    if (this->dict().isDict("powerOffCriterionModel"))
    {
        powerOffCriterionModelPtr_.reset
        (
            powerOffCriterionModel::New
            (
                mesh,
                this->dict().subDict("powerOffCriterionModel")
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::structureModel::correct
(
    const volScalarField& HT,
    const volScalarField& H
)
{
    //- Correct powerModels & correct active structure surface temperature
    forAllIter
    (
        powerModelTable,
        powerModels_,
        iter
    )
    {
        iter()->correct(HT, H);
        iter()->correctT(Tact_);
    }
    Tact_.correctBoundaryConditions();

    if (Tpas_.writeOpt() == IOobject::AUTO_WRITE)
    {
        //- Correct inert subStructure. What follows is the equivalent of doing
        //  the following:
        /*
            fvScalarMatrix pasEqn
            (
                    fvm::ddt(alphaRhoCppas_, Tpas_)
                ==
                    iApas_*HT
                -   fvm::Sp(iApas_*H, Tpas_)
            );
            pasEqn.solve();
        */
        //  Except, it is faster like this rather than to solve an equation 
        //  over the entire mesh, as the passive subStructure might not exist
        //  everywhere

        scalar dt(mesh_.time().deltaT().value());
        const volScalarField& Tpas0(Tpas_.oldTime()); 
        forAll(cells_, i)
        {
            label celli(cells_[i]);
            const scalar& iA(iApas_[celli]);
            if (iA == 0) continue;  //- Avoid solving where the passive 
                                    //  structure does not exist
            scalar alphaRhoCpByDt(alphaRhoCppas_[celli]/dt);
            Tpas_[celli] = 
                (
                    iA*HT[celli] 
                +   alphaRhoCpByDt*Tpas0[celli]
                )/
                (alphaRhoCpByDt + iA*H[celli]);
        }
        Tpas_.correctBoundaryConditions();
    }

    //- Update Twall, heatFlux    
    forAll(cells_, i)
    {
        const label& celli(cells_[i]);
        scalar& Twall(Twall_[celli]);

        //- Set indicative wall temperature as max between power structure
        //  surface temperature and passive structure surface temperature
        Twall = Foam::max(Tact_[celli], Tpas_[celli]);

        //- Update heat flux (mostly for extra info purposes, maybe only
        //  used by the Shah pool boiling model under some circumstances).
        //  For representativity, it is set as the max heat flux to the fluid
        //  calculated between the active and passive subStructures
        heatFlux_[celli] = H[celli]*Twall-HT[celli];
    }
    Twall_.correctBoundaryConditions();
    heatFlux_.correctBoundaryConditions();
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

void Foam::structureModel::checkPowerOff()
{
    if (powerOffCriterionModelPtr_.valid())
    {
        if (powerOffCriterionModelPtr_->powerOffCriterion())
        {
            forAllIter
            (
                powerModelTable,
                powerModels_,
                iter
            )
            {
                iter()->powerOff();
            }
        }
    }
}

// ************************************************************************* //
