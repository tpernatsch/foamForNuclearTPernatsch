/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "twoPhase.H"
#include "addToRunTimeSelectionTable.H"
#include "myMULES.H"
#include "regime.H"
#include "regimeMapModel.H"
#include "saturationModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalHydraulicModels
{
    defineTypeNameAndDebug(twoPhase, 0);
    addToRunTimeSelectionTable
    (
        thermalHydraulicModel, 
        twoPhase, 
        thermalHydraulicModels
    );

    const Enum<twoPhase::phaseChangeModel>
    phaseChangeModelNames_
    (
        {
            { 
                twoPhase::phaseChangeModel::none, 
                "none" 
            },
            { 
                twoPhase::phaseChangeModel::heatConductionLimited, 
                "heatConductionLimited" 
            },
            { 
                twoPhase::phaseChangeModel::individualSideDriven, 
                "individualSideDriven" 
            },
        }
    );
}
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalHydraulicModels::twoPhase::twoPhase
(
    const fvMesh& mesh,
    pimpleControl& pimple,
    fv::options& fvOptions
)
:
    thermalHydraulicModel
    (
        mesh,
        pimple,
        fvOptions
    ),
    fluid1_
    (
        this->subDict
        (
                word(this->lookup("fluid1"))
            +   "Properties"
        ),
        mesh,
        word(this->lookup("fluid1"))
    ),
    fluid2_
    (
        this->subDict
        (
                word(this->lookup("fluid2"))
            +   "Properties"
        ),
        mesh,
        word(this->lookup("fluid2"))
    ),
    structurePtr_
    (
        structureModel::New
        (
            this->subDict("structureProperties"),
            mesh   
        )
    ),
    structure_(structurePtr_()),
    movingAlpha_
    (
        1.0 - structure_
    ),
    FFPair_(fluid1_, fluid2_),
    F1SPair_(fluid1_, structure_),
    F2SPair_(fluid2_, structure_),
    residualKd_
    (
        dimensionedScalar::lookupOrDefault
        (
            "residualKd",
            *this,
            dimMass/dimVol/dimTime,
            1
        )
    ),
    iA12_
    (
        IOobject
        (
            "iA."+fluid1_.name()+"."+fluid2_.name(),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimArea/dimVol, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    frac1_
    (
        IOobject
        (
            "frac."+fluid1_.name(),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    frac2_
    (
        IOobject
        (
            "frac."+fluid1_.name(),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    dmdt_
    (
        IOobject
        (
            "dmdt",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimMass/dimVol/dimTime, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    iT12_
    (
        IOobject
        (
            "T.interface",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    regimeMap_
    (
        regimeMapModel::New
        (
            mesh,
            this->subDict("regimeMapModel"),
            this->subDict("physicsModelsByRegime")
        )
    ),
    phaseChangeModel_
    (
        phaseChangeModelNames_.get
        (
            this->lookupOrDefault<word>("phaseChangeModel", "none")
        )
    ),
    phaseChange_
    (
        (phaseChangeModel_ != phaseChangeModel::none) ?
        true :
        false
    ),
    volatilePhaseName_
    (
        this->lookupOrDefault<word>("volatilePhase", fluid2_.name())
    ),
    partialElimination_
    (
        pimple.dict().lookupOrDefault<Switch>("partialElimination", false)
    ),
    faceMomentum_
    (
        pimple.dict().lookupOrDefault<Switch>("faceMomentum", false)
    )
{
    //- Normalize phase fraction fields, structure has priority over fluid2, 
    //  fluid2 has priority over fluid1
    fluid1_.volScalarField::operator=(1.0-fluid2_);
    volScalarField corr(movingAlpha_/(fluid1_+fluid2_));
    fluid1_.volScalarField::operator*=(corr);
    fluid2_.volScalarField::operator*=(corr);

    //- Initialize phi
    phi_ = 
        fvc::interpolate(fluid1_)*fluid1_.phi()
      + fvc::interpolate(fluid2_)*fluid2_.phi();
    
    //- Initialize dragCoefficient and heatTransferCoefficient tables
    #include "initTables.H"

    if (phaseChange_)
    {
        saturation_.reset
        (
            saturationModel::New
            (
                this->subDict("saturationModel"),
                fluid1_,
                fluid2_
            )
        );
    }

    //- Set initial interfacial temperature (to avoid problems with the
    //  under-relaxation at the first time-step in EEqns.H
    IOobject iT12Header
    (
        iT12_.name(),
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ
    );
    if (!iT12Header.typeHeaderOk<volScalarField>(true))
    {
        if (phaseChange_)
        {
            iT12_ = saturation_->Tsat(p_);
        }
    }

    if (phaseChangeModel_ == phaseChangeModel::individualSideDriven)
    {
        if  (
                volatilePhaseName_ != fluid1_.name() and 
                volatilePhaseName_ != fluid2_.name()
            )
        {
            FatalErrorInFunction
            << "Phase " << volatilePhaseName_ << " not found!"
            << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::thermalHydraulicModels::twoPhase::rho() const
{
    return fluid1_*fluid1_.thermo().rho() + fluid2_*fluid2_.thermo().rho();
}

Foam::tmp<Foam::volVectorField> 
Foam::thermalHydraulicModels::twoPhase::U() const
{
    return fluid1_*fluid1_.U() + fluid2_*fluid2_.U();
}

void Foam::thermalHydraulicModels::twoPhase::correct(scalar& residual)
{   
    #include "solve.H"
    
    Info << endl;
}

void Foam::thermalHydraulicModels::twoPhase::correctCourant()
{
    CoNum_ = 0.0;
    meanCoNum_ = 0.0;

    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi_))().primitiveField()
    );

    CoNum_ = 0.5*gMax(sumPhi/mesh_.V().field())*runTime_.deltaTValue();

    meanCoNum_ =
        0.5*(gSum(sumPhi)/gSum(mesh_.V().field()))*runTime_.deltaTValue();

    Info<< "Courant Number (avg max) = " << meanCoNum_
        << " " << CoNum_ << endl;

    scalar UrCoNum = 0.5*gMax
    (
            fvc::surfaceSum
            (
                    mag(fluid1_.phi() 
                -   fluid2_.phi())
            )().primitiveField()
        /   mesh_.V().field()
    )*runTime_.deltaTValue();

    Info<< "Max Ur Courant Number = " << UrCoNum << endl;

    CoNum_ = max(CoNum_, UrCoNum);
}

void Foam::thermalHydraulicModels::twoPhase::correctContErrs()
{
    volScalarField& rho1(fluid1_.thermo().rho());
    volScalarField& rho2(fluid2_.thermo().rho());
    
    fluid1_.contErr() = 
    (
            fvc::ddt(fluid1_, rho1)
        +   fvc::div(fluid1_.alphaRhoPhi())
        -   (fvOptions_(fluid1_, rho1) & rho1)
        + dmdt_ 
    );
    fluid1_.contErr().correctBoundaryConditions();

    fluid2_.contErr() = 
    (
            fvc::ddt(fluid2_, rho2) 
        +   fvc::div(fluid2_.alphaRhoPhi())
        -   (fvOptions_(fluid2_, rho2) & rho2)
        - dmdt_
    );
    fluid2_.contErr().correctBoundaryConditions();
}


// ************************************************************************* //