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
#include "phaseChangeModel.H"
#include "myOps.H"
//#include <chrono>

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
}
}

const Foam::Enum
<
    Foam::thermalHydraulicModels::twoPhase::partialEliminationMode
>
Foam::thermalHydraulicModels::twoPhase::partialEliminationModeNames_
(
    {
        { 
            partialEliminationMode::none, 
            "none" 
        },
        { 
            partialEliminationMode::legacy, 
            "explicit" 
        },
        { 
            partialEliminationMode::implicit, 
            "implicit" 
        }
    }
);

const Foam::Enum
<
    Foam::thermalHydraulicModels::twoPhase::contErrCompensationMode
>
Foam::thermalHydraulicModels::twoPhase::contErrCompensationModeNames_
(
    {
        { 
            contErrCompensationMode::none, 
            "none" 
        },
        { 
            contErrCompensationMode::Su, 
            "Su" 
        },
        { 
            contErrCompensationMode::Sp, 
            "Sp" 
        },
        { 
            contErrCompensationMode::SuSp, 
            "SuSp" 
        }
    }
);

const Foam::Enum
<
    Foam::thermalHydraulicModels::twoPhase::heStabilizationMode
>
Foam::thermalHydraulicModels::twoPhase::heStabilizationModeNames_
(
    {
        { 
            heStabilizationMode::cutoff, 
            "cutoff" 
        },
        { 
            heStabilizationMode::source, 
            "source" 
        }
    }
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class fieldType>
void Foam::thermalHydraulicModels::twoPhase::relaxFieldPtrTable
(
    const word& defaultRelaxFactorName,
    HashPtrTable
    <
        fieldType,
        word, 
        word::hash
    >& fieldPtrTable
)
{
    scalar f0(myOps::relaxationFactor(mesh_, defaultRelaxFactorName));    
    forAll(fieldPtrTable.toc(), i)
    {   
        word key(fieldPtrTable.toc()[i]);
        fieldType& field(*fieldPtrTable[key]);
        //- If a relaxation factor was provided specifically for one of the
        //  table fields, use that instead of the one specified by 
        //  defaultRelaxFactorName
        scalar f = 
            (myOps::relax(mesh_, field.name())) ?
            myOps::relaxationFactor(mesh_, field.name()) :
            f0;

        if (f != 1.0)
        {
            field = f*field+(1.0-f)*field.prevIter();
            field.correctBoundaryConditions();
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalHydraulicModels::twoPhase::twoPhase
(
    Time& time,
    fvMesh& mesh,
    myPimpleControl& pimple,
    fv::options& fvOptions
)
:
    thermalHydraulicModel
    (
        time,
        mesh,
        pimple,
        fvOptions
    ),
    //- The structure model needs to be created before the fluids because
    //  the turbulence models created by the fluids might require a reference
    //  to the structure (obtained via objectRegistry lookup in the specific
    //  turbulence model class)
    structurePtr_
    (
        structureModel::New
        (
            this->subDict("structureProperties"),
            mesh   
        )
    ),
    structure_(structurePtr_()),
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
    movingAlpha_
    (
        "movingAlpha",
        1.0 - structure_
    ),
    FFPair_(fluid1_, fluid2_),
    F1SPair_(fluid1_, structure_),
    F2SPair_(fluid2_, structure_),
    ReTwoPhase_
    (
        IOobject
        (
            "Re.mixture.structure",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 10)
    ),
    U_
    (
        IOobject
        (
            "U",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("", dimVelocity, vector::zero)
    ),
    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimDensity, 0.0)
    ),
    Vm_(this->lookupOrDefault<scalar>("virtualMassCoeff", 0)),
    iA12_
    (
        IOobject
        (
            "areaDensity.interface",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimArea/dimVol, 0),
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
    VmVolForce1_
    (
        fluid1_.U(), dimMass*dimLength/dimTime/dimTime
    ),
    VmVolForce2_
    (
        fluid2_.U(), dimMass*dimLength/dimTime/dimTime
    ),
    cumulContErr1_(0),
    cumulContErr2_(0),
    bothPhasesArePresent_(false),
    withinMarginToPhaseChange_(false),
    partialEliminationMode_
    (
        partialEliminationModeNames_.get
        (
            pimple_.dict().lookupOrDefault<word>
            (
                "partialEliminationMode", 
                "none"
            )
        )
    ),
    UContErrCompensationMode_
    (
        contErrCompensationModeNames_.get
        (
            (pimple_.dict().found("UContinuityErrorCompensationMode")) ?
                pimple_.dict().get<word>
                (
                    "UContinuityErrorCompensationMode"
                ) :
                pimple_.dict().lookupOrDefault<word>
                (
                    "continuityErrorCompensationMode",
                    "SuSp"
                )
        )
    ),
    heContErrCompensationMode_
    (
        contErrCompensationModeNames_.get
        (
            (pimple_.dict().found("heContinuityErrorCompensationMode")) ?
                pimple_.dict().get<word>
                (
                    "heContinuityErrorCompensationMode"
                ) :
                pimple_.dict().lookupOrDefault<word>
                (
                    "continuityErrorCompensationMode",
                    "SuSp"
                )
        )
    ),
    heStabilizationMode_
    (
        heStabilizationModeNames_.get
        (
            pimple_.dict().lookupOrDefault<word>
            (
                "enthalpyStabilizationMode",
                "cutoff"
            )
        )
    ),
    oscillationLimiterFraction_
    (
        pimple.dict().lookupOrDefault<scalar>
        (
            "oscillationLimiterFraction", 
            0.0
        )
    )
{
    //- Init autoPtr-managed fields, namely flowQuality, XLM (i.e. 
    //  Lockhart-Martinelli parameter) and dispersion
    fluid1_.initTwoPhaseFields();
    fluid2_.initTwoPhaseFields();
    
    //- This is specifically to avoid problem when solveAlpha is set to solve
    //  for only one phase, the one for which BCs and initial conditions are
    //  provided, yet this phase starts at 0 while the other phase BC and
    //  initial conditions were not provided. This is problematic as by 
    //  default a phase defaults to 0, so if the provided phase is 0 everywhere
    //  too, you have issues with the normalization later. This is veeeery
    //  specific but hey, why not!
    IOobject fluid1Header
    (
        "alpha."+fluid1_.name(),
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ
    );
    IOobject fluid2Header
    (
        "alpha."+fluid2_.name(),
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ
    );
    if 
    (
        !fluid1Header.typeHeaderOk<volScalarField>(true)
    and fluid2Header.typeHeaderOk<volScalarField>(true)
    )
    {
        fluid1_.volScalarField::operator=(geometricOneField()-fluid2_);
    }
    if 
    (
        !fluid2Header.typeHeaderOk<volScalarField>(true)
    and fluid1Header.typeHeaderOk<volScalarField>(true)
    )
    {
        fluid2_.volScalarField::operator=(geometricOneField()-fluid1_);
    }

    //- Normalize phase fraction fields, structure is left unchanged
    volScalarField corr(movingAlpha_/(fluid1_+fluid2_));
    fluid1_.volScalarField::operator*=(corr);
    fluid2_.volScalarField::operator*=(corr);
    fluid1_.correctBoundaryConditions();
    fluid2_.correctBoundaryConditions();

    //- Set the normalized phase fraction fields
    fluid1_.normalized() = fluid1_/movingAlpha_;
    fluid2_.normalized() = fluid2_/movingAlpha_;

    //- Set the bothPhasesArePresentFlag. This is only used within the
    //  adjustTimeStep function and for avoiding doing subcycles if
    //  only one phase is present
    bothPhasesArePresent_ = 
    (
        max((fluid1_*fluid2_)()).value() >= 1e-4 
    );

    //- Initialize fluid-intensive fluxes (i.e. that depend on the phase
    //  fraction, namely alphaPhi and alphaRhoPhi, which are the REAL 
    //  volumetric flux in m3/s and the REAL mass flux in kg/s. By REAL I mean
    //  not superficial). This is done after the phaseFraction normalization 
    //  step to ensure consistency. This step has an effect ONLY IF the
    //  alphaPhi, alphaRhoPhi fields were NOT found on disk
    fluid1_.initAlphaPhis();
    fluid2_.initAlphaPhis();

    //- Set total volumetric flux (real one, not superficial)
    phi_ = fluid1_.alphaPhi() + fluid2_.alphaPhi();

    //- Initialize drag coeff and heat transfer coeff tables
    #include "initTables_2p.H"

    //- Create turbulence models. This is done outside of fluid constructors as
    //  turbulence models might require references to fields that do not exist
    //  yet (e.g. the Reynolds number between fluid and structure, or between
    //  the two fluids, which are in the FFPair/FSPair classes, that cannot be
    //  created before the fluids)
    fluid1_.constructTurbulenceModel();
    fluid2_.constructTurbulenceModel();

    //- Initialize local drag coeff and heat transfer coeff tables for each
    //  (non-interpolated) regime
    forAllIter
    (
        regimeTable,
        regimeMap_->regimes(),
        regimeIter
    )
    {
        regime& regime(regimeIter()());
        regime.initLocalTables(Kds_, htcs_);
    }

    //- Init phaseChangeModel
    if (this->found("phaseChangeModel"))
    {
        phaseChange_.reset
        (
            phaseChangeModel::New
            (
                this->subDict("phaseChangeModel"),
                pimple_,
                fluid1_,
                fluid2_,
                p_,
                htcs_,
                dmdt_,
                iT12_,
                iA12_
            )
        );
    }

    //-
    U_ = this->U();
    rho_ = this->rho();

    //- Compute initialFluidMass (only relevant for compressibility effects,
    //  yet it is still work in progress)
    initialFluidMass_ = fvc::domainIntegrate(this->rho()*movingAlpha_);

    Info << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::thermalHydraulicModels::twoPhase::rho() const
{
    return 
        (
            fluid1_*fluid1_.thermo().rho() 
        +   fluid2_*fluid2_.thermo().rho()
        )/movingAlpha_;
}


Foam::tmp<Foam::volVectorField> 
Foam::thermalHydraulicModels::twoPhase::U() const
{
    return fluid1_*fluid1_.U() + fluid2_*fluid2_.U();
}


void Foam::thermalHydraulicModels::twoPhase::correct
(
    scalar& residual,
    bool solveFluidDynamics, 
    bool solveEnergy
)
{   
    correctRegimes(solveFluidDynamics, solveEnergy);
    if (solveFluidDynamics)
    {
        correctFluidMechanics(residual);
    }
    
    if (solveEnergy)
    {
        correctEnergy(residual);
    }

    Info << endl;
}


void Foam::thermalHydraulicModels::twoPhase::correctEnergy(scalar& residual)
{
    //- Solve fluid energy equations and update structure energy sources
    #include "EEqns_2p.H"
}


void Foam::thermalHydraulicModels::twoPhase::correctFluidMechanics
(
    scalar& residual
)
{
    correctContErrs();

    //- Solve continuity equations to compute new alphas
    #include "alphaEqns_2p.H"
    
    //- Construct momentum matrices
    #include "UEqns_2p.H"

    //- Solve pressure equation and reconstruct velocities
    if (momentumMode_ == momentumMode::faceCentered)
    {
        #include "pEqnf_2p.H"
    }
    else
    {
        #include "pEqn_2p.H"
    }
    
    //- Continuity error adjustments and infos
    correctContErrs();
    printContErrs();
    calcCumulContErrs();
}


void Foam::thermalHydraulicModels::twoPhase::correctRegimes
(
    bool solveFluidDynamics, 
    bool solveEnergy
)
{
    //- Update spatial extent of regimes
    regimeMap_->correct();
    
    //- Reset all relevant fields that need to be set by the regimes
    forAllIter
    (
        volTensorFieldPtrTable,
        Kds_,
        iter
    )
    {
        volTensorField& Kd(*iter());
        //- Cache previous field for underrelaxation
        Kd.storePrevIter(); 
        Kd *= 0.0;
    }
    forAllIter
    (
        volScalarFieldPtrTable,
        htcs_,
        iter
    )
    {
        volScalarField& htc(*iter());
        //- Cache previous field of underrelaxation
        //*oldHtcs_[iter.key()] = htc;
        htc.storePrevIter();
        htc *= 0.0;
    }

    //- Update continuity errors, incluisve of mass transfer. As mass transfer
    //  can vary also if not doing fluidMechanics, continuityErrors are always
    //  corrected, just in case (even though I have no clue why one would
    //  do phase change simulations with fluidMechanics off, maybe a weird 
    //  steady state?)
    //correctContErrs();

    //- Correct local regime fluid geometry models
    forAllIter
    (
        regimeTable,
        regimeMap_->regimes(),
        regimeIter
    )
    {
        regime& regime(regimeIter()());
        if (!regime.requiresModelCorrection()) continue;
        regime.correctFluidGeometryModels();
    }

    //- Some fluidGeometry models (e.g. sphericalTopology of 
    //  fluidInterfacialArea) might depend on other fluidGeometry quantities,
    //  so reset them only after having updated the local regime fluidModel
    fluid1_.dispersion() *= 0.0;
    fluid2_.dispersion() *= 0.0;
    fluid1_.Dh() *= 0.0;
    fluid2_.Dh() *= 0.0;

    if(mesh_.relaxField(iA12_.name()))
        iA12_.storePrevIter();
    iA12_ *= 0.0;

    //- Correct global fluid geometry fields. Why didn't I do these two loops
    //  in one? Well, in the loop below, if the regime is interpolated, It will
    //  be looking at the respective fields of interest (iA, fracs, etc.) of 
    //  both the regimes from which it is interpolated from. However, nothing
    //  guarantees that by the time I'm correcting an interpolated regime, I am
    //  done correcting the regimes from which it is interpolated from. Thus, 
    //  the previous loop corrects all non-interpolated regimes before using
    //  them (among the non-interpolated ones) to correct the global fields
    forAllConstIter
    (
        regimeTable,
        regimeMap_->regimes(),
        regimeIter
    )
    {
        const regime& regime(regimeIter());
        if (!regime.isCurrentlyPresent()) continue;
        regime.correctFluidGeometryFields
        (
            iA12_,
            fluid1_,
            fluid2_
        );
    }

    //- Correct BCs of global fields outside of regimes, as not the same
    //  regimes might be present on all processors, meaning that a 
    //  correctBoundaryConditions would be called on a global field an unequal
    //  amount of times, resulting in an MPI wait error. Nonetheless, the
    //  usefuleness of these BCs is debatable given that these fields are
    //  (amlost?) always maniupulated on a cell-by-cell basis, except maybe for
    //  the Dh
    iA12_.correctBoundaryConditions();
    fluid1_.Dh().correctBoundaryConditions();
    fluid2_.Dh().correctBoundaryConditions();

    //- Relax iA12_
    iA12_.relax();

    //- Update dimensionless numbers (Re, etc.). Done now as the previous loop
    //  corrected (among others) the fluid characteristic dimensions which are
    //  needed e.g. to compute the Reynolds numbers
    FFPair_.correct();
    F1SPair_.correct();
    F2SPair_.correct();
    ReTwoPhase_ = 
        max
        (
            mag
            (
                fluid1_.normalized()*fluid1_.rho()*fluid1_.U()
            +   fluid2_.normalized()*fluid2_.rho()*fluid2_.U()
            )*structure_.Dh()/
            (
                fluid1_.normalized()*fluid1_.thermo().mu()
            +   fluid2_.normalized()*fluid2_.thermo().mu()
            ),
            dimensionedScalar("", dimless, 10)
        );       
    //- Correct local drag and heat transfer models. As these might depend on
    //  dimensionless numbers such as Re, this must be done only after having
    //  updated the pairs, which in turn is done only after updating the fluid
    //  geometry. This is why it cannot be done in one loop
    forAllIter
    (
        regimeTable,
        regimeMap_->regimes(),
        regimeIter
    )
    {
        regime& regime(regimeIter()());
        if (!regime.requiresModelCorrection()) continue;
        if (solveFluidDynamics) regime.correctDragModels();
        if (solveEnergy) regime.correctHeatTransferModels();
    }

    //- Correct global drag and heat transfer coefficients' table. As stated
    //  before, this loop and the previous one cannot be merged as nothing
    //  guarantees that the local drag/heat transfer models of all the non-
    //  interpolated regimes will have been updated before updating the
    //  local interpolated regimes
    forAllConstIter
    (
        regimeTable,
        regimeMap_->regimes(),
        regimeIter
    )
    {
        const regime& regime(regimeIter());
        if (!regime.isCurrentlyPresent()) continue;
        if (solveFluidDynamics) regime.correctDragTable(Kds_);
        if (solveEnergy) regime.correctHeatTransferTable(htcs_); 
    }

    //- Correct BCs of global fields outside of regimes, as not the same
    //  regimes might be present on all processors, meaning that a 
    //  correctBoundaryConditions would be called on a global field an unequal
    //  amount of times, resulting in an MPI wait error. Nonetheless, the
    //  usefuleness of these BCs is debatable given that these fields are
    //  (amlost?) always maniupulated on a cell-by-cell basis
    forAllIter
    (
        volTensorFieldPtrTable,
        Kds_,
        iter
    )
    {
        volTensorField& Kd(*iter());
        Kd.correctBoundaryConditions();
    }
    forAllIter
    (
        volScalarFieldPtrTable,
        htcs_,
        iter
    )
    {
        volScalarField& htc(*iter());
        htc.correctBoundaryConditions();
    }
}


void Foam::thermalHydraulicModels::twoPhase::correctCourant()
{
    CoNum_ = 0.0;
    meanCoNum_ = 0.0;

    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi_))().primitiveField()/
        movingAlpha_.primitiveField()
    );

    CoNum_ = 0.5*gMax(sumPhi/mesh_.V().field())*runTime_.deltaTValue();

    meanCoNum_ =
        0.5*(gSum(sumPhi)/gSum(mesh_.V().field()))*runTime_.deltaTValue();

    Info<< "Courant Number (avg max) = " << meanCoNum_
        << " " << CoNum_ << endl;

    scalar UrCoNum = 
        0.5*gMax
        (
            fvc::surfaceSum
            (
                mag(fluid1_.phi()-fluid2_.phi())
            )().primitiveField()/
            movingAlpha_/
            mesh_.V().field()
        )*runTime_.deltaTValue();

    Info<< "Max Ur Courant Number = " << UrCoNum << endl;

    CoNum_ = max(CoNum_, UrCoNum);
}


void Foam::thermalHydraulicModels::twoPhase::adjustTimeStep()
{
    this->correctCourant();

    bool phase1(max(fluid1_).value() >= 1e-6);
    bool phase2(max(fluid2_).value() >= 1e-6);
    bothPhasesArePresent_ = (phase1 and phase2);

    bool adjustTimeStep =
        runTime_.controlDict().lookupOrDefault("adjustTimeStep", false);

    if (adjustTimeStep)
    {
        scalar maxCo =
            runTime_.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

        if (runTime_.controlDict().found("maxCoTwoPhase"))
        {
            scalar maxCoTwoPhase
            (
                runTime_.controlDict().get<scalar>("maxCoTwoPhase")
            );

            if (bothPhasesArePresent_) maxCo = maxCoTwoPhase;
            else
            {
                if 
                (
                    runTime_.controlDict().found("marginToPhaseChange")
                and phaseChange_.valid()
                )
                {
                    scalar marginToPhaseChange
                    (
                        runTime_.controlDict().get<scalar>
                        (
                            "marginToPhaseChange"
                        )
                    );

                    scalar DT1(min(mag(fluid1_.T()-iT12_)().primitiveField()));
                    scalar DT2(min(mag(fluid2_.T()-iT12_)().primitiveField()));
                    if 
                    (
                        (
                            phase1 and !phase2 and DT1 < marginToPhaseChange
                        ) or
                        (
                            phase2 and !phase1 and DT2 < marginToPhaseChange
                        )
                    )
                    {
                        maxCo = maxCoTwoPhase;
                    }
                }
            }
        }

        scalar maxDeltaT =
            runTime_.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

        scalar minDeltaT =
            runTime_.controlDict().lookupOrDefault<scalar>("minDeltaT", 1e-69);

        scalar f =
            std::abs
            (
                runTime_.controlDict().lookupOrDefault<scalar>
                (
                    "maxDeltaTMaxRelInc", 
                    0.1
                )
            );
        scalar maxDeltaTFact = maxCo/(CoNum_ + SMALL);
        scalar deltaTFact = 
            min
            (
                min
                (
                    maxDeltaTFact, 
                    1.0 + f*maxDeltaTFact
                ), 
                1.0 + f
            );

        runTime_.setDeltaT
        (
            max
            (
                min
                (
                    deltaTFact*runTime_.deltaTValue(),
                    maxDeltaT
                ),
                minDeltaT
            )
        );
    }
}


void Foam::thermalHydraulicModels::twoPhase::correctContErrs()
{
    volScalarField& cE1(fluid1_.contErr());
    volScalarField& cE2(fluid2_.contErr());
    volScalarField& rho1(fluid1_.thermo().rho());
    volScalarField& rho2(fluid2_.thermo().rho());

    cE1 = 
    (
        fvc::ddt(fluid1_, rho1)
    +   fvc::div(fluid1_.alphaRhoPhi())
    -   (fvOptions_(fluid1_, rho1) & rho1)
    +   dmdt_
    );
    cE2 = 
    (
        fvc::ddt(fluid2_, rho2) 
    +   fvc::div(fluid2_.alphaRhoPhi())
    -   (fvOptions_(fluid2_, rho2) & rho2)
    -   dmdt_
    );
    
    cE1.correctBoundaryConditions();
    cE2.correctBoundaryConditions();
}


void Foam::thermalHydraulicModels::twoPhase::printContErrs()
{
    volScalarField contErrRel1(fluid1_.contErr()/fluid1_.rho());
    volScalarField contErrRel2(fluid2_.contErr()/fluid2_.rho());

    Info<< "Instantaneous relative continuity errors (" 
        << fluid1_.name() << " " << fluid2_.name() << ") = "
        << contErrRel1.weightedAverage(mesh_.V()).value() << " "
        << contErrRel2.weightedAverage(mesh_.V()).value() << " "
        << "1/s" << endl;
}

void Foam::thermalHydraulicModels::twoPhase::calcCumulContErrs()
{
    if (pimple_.finalIter())
    {
        scalarField integralContErr1
        (
            mesh_.time().deltaTValue()*
            fvc::volumeIntegrate(fluid1_.contErr())
        );

        scalarField integralContErr2
        (
            mesh_.time().deltaTValue()*
            fvc::volumeIntegrate(fluid2_.contErr())
        );

        const scalarField& V(mesh_.V());
        scalar totV(0);

        forAll(V, i)
        {
            cumulContErr1_ += integralContErr1[i];
            cumulContErr2_ += integralContErr2[i];
            totV += V[i];
        }

        Info<< "Cumulative continuity errors ("
            << fluid1_.name() << " " << fluid2_.name() << ") = " 
            << (cumulContErr1_/totV) << " " << (cumulContErr2_/totV)
            << " kg/m3" << endl;
    }
}

// ************************************************************************* //