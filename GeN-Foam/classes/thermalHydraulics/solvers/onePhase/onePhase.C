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

#include "onePhase.H"
#include "addToRunTimeSelectionTable.H"
#include "regime.H"
#include "regimeMapModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalHydraulicModels
{
    defineTypeNameAndDebug(onePhase, 0);
    addToRunTimeSelectionTable
    (
        thermalHydraulicModel, 
        onePhase, 
        thermalHydraulicModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalHydraulicModels::onePhase::onePhase
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
    fluid_
    (
        this,
        mesh,
        word(""),   //- This is the phase name, setting it to "" signals a
                    //  onePhase solver to the rest of the FFS library
        false       //- No need to read or write the fluid phaseFraction in
                    //  onePhase, it is tied to the structure phaseFraction
    ),
    FSPair_(fluid_, structure_),
    regimeMap_
    (
        regimeMapModel::New
        (
            mesh,
            this->subDict("regimeMapModel"),
            this->subDict("physicsModelsByRegime")
        )
    )
{
    //- Create turbulence model
    fluid_.constructTurbulenceModel();

    //- Normalize phase fraction fields, structure has priority
    fluid_.volScalarField::operator=(1.0-structure_);

    //- Set fluid characteristic dimension to structure hydraulic diameter.
    //  This is handled by the fluidGeometry class in the twoPhase solver
    //  (as the fluid characteristic dimension will depend on the regime)
    //  so here I need to do it manually. Since I assume that the structure is
    //  immutable, this is done only once
    fluid_.Dh() = structure_.Dh();

    //- Initialize fluid-intensive fluxes (i.e. that depend on the phase
    //  fraction, namely alphaPhi and alphaRhoPhi, which are the REAL 
    //  volumetric flux in m3/s and the REAL mass flux in kg/s. By REAL I mean
    //  not superficial). This is done after the phaseFraction normalization 
    //  step to ensure consistency. This step has an effect ONLY IF the
    //  alphaPhi, alphaRhoPhi fields were NOT found on disk
    fluid_.initAlphaPhis();

    //- Initialize continuity errors. It's important to do it here or, if not
    //  solving for fluidMechanics, these would never get corrected
    correctContErr();

    //- The total volumetric flux is the REAL fluid volumetric flux. Might
    //  as well remove phi_ entirely as a field, I know... Maybe in the future
    phi_ = fluid_.alphaPhi();

    //- Initialize dragCoefficient and heatTransferCoefficient tables.
    //  These tables actually consist of only one entry coresponding to the
    //  fluid-structure pair. Why have a table then?  Well, I first
    //  developed the twoPhase solver, which requires tables as there are
    //  multiple pairs. Since I want to keep the same base library and code
    //  structure for both the single and two phase, I decided to keep this
    //  rather than unnecessarily duplicate the base code of the models library
    Kds_.insert
    (
        "FSPair",
        autoPtr<volTensorField>
        (
            new volTensorField
            (
                IOobject
                (
                    "Kd",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedTensor
                (
                    "", 
                    dimDensity/dimTime, 
                    tensor(0,0,0,0,0,0,0,0,0)
                ),
                zeroGradientFvPatchScalarField::typeName
            )
        )
    );

    htcs_.insert
    (
        "FSPair",
        autoPtr<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "htc",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("", dimPower/dimArea/dimTemperature, 0),
                zeroGradientFvPatchScalarField::typeName
            )
        )
    );

    //- Initialize local drag coeff and heat transfer coeff tables for each
    //  (non-interpolated) regime. 
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

    //- Compute initialFluidMass
    initialFluidMass_ = fvc::domainIntegrate(fluid_.rho()*fluid_);

    Info << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- Solve according to flags
void Foam::thermalHydraulicModels::onePhase::correct
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

void Foam::thermalHydraulicModels::onePhase::correctFluidMechanics
(
    scalar& residual
)
{
    #include "UEqn_1p.H"
    if (momentumMode_ == momentumMode::faceCentered)
    {
        #include "pEqnf_1p.H"
    }
    else
    {
        #include "pEqn_1p.H"
    }

    //- Continuity error adjustment and infos
    correctContErr();
    printContErr();
    calcCumulContErr();
}

void Foam::thermalHydraulicModels::onePhase::correctEnergy(scalar& residual)
{
    #include "EEqn_1p.H"
}

void Foam::thermalHydraulicModels::onePhase::correctRegimes
(
    bool solveFluidDynamics, 
    bool solveEnergy
)
{
    //- Reset all fields relevant for fluid-structure coupling
    forAllIter
    (
        volTensorFieldPtrTable,
        Kds_,
        iter
    )
    {
        volTensorField& Kd(*iter());
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
        htc *= 0.0;
    }

    //- Update continuity errors. These only depend on alphaPhi and alphaRhoPhi
    //  so its useless to updated them in not doing solveFluidDynamics, as they
    //  would be constant in such case
    if (solveFluidDynamics) correctContErr();

    //- Correct regimes (i.e. regime marker fields)
    regimeMap_->correct();

    //- Update dimensionless numbers (Re, etc.)
    FSPair_.correct();

    //- Correct field tables for drag and heat transfer according
    //  to new regime distribution in domain
    //- It is kind of overkill as the table consist only of one entry,
    //  but I do not want to duplicate the code by creating single-phase
    //  variants of the correctDragTable and correctHeatTransferTable
    //  functions, that'd be pointless and the overhead associated with
    //  looping over a hashTable with a single entry is irrelevant

    //- Correct local drag and heat transfer models. As these might depend on
    //  dimensionless numbers such as Re, this must be done only after having
    //  updated the FSPair that contains the dimensionless numbers that might
    //  be used by drag or heatTransfer models
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

    //- Correct global drag and heat transfer coefficients' table. This loop
    //  and the previous one cannot be merged as nothing guarantees that the
    //  local drag/heat transfer models of all the non-interpolated regimes
    //  will have been updated before updating the local interpolated regimes
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
}

void Foam::thermalHydraulicModels::onePhase::correctCourant()
{
    CoNum_ = 0.0;
    meanCoNum_ = 0.0;

    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi_))().primitiveField()/fluid_.primitiveField()
    );

    CoNum_ = 0.5*gMax(sumPhi/mesh_.V().field())*runTime_.deltaTValue();

    meanCoNum_ =
        0.5*(gSum(sumPhi)/gSum(mesh_.V().field()))*runTime_.deltaTValue();

    Info<< "Courant Number (avg max) = " << meanCoNum_
        << " " << CoNum_ << endl;
}

void Foam::thermalHydraulicModels::onePhase::correctContErr()
{
    volScalarField& cE(fluid_.contErr());
    volScalarField& rho(fluid_.rho());

    cE = 
    (
        fvc::ddt(fluid_, rho)
    +   fvc::div(fluid_.alphaRhoPhi())
    -   (fvOptions_(fluid_, rho) & rho)
    );
    
    cE.correctBoundaryConditions();
}


void Foam::thermalHydraulicModels::onePhase::printContErr()
{
    volScalarField contErrRel(fluid_.contErr()/fluid_.rho());

    Info<< "Instantaneous relative continuity error (avg) = "
        << contErrRel.weightedAverage(mesh_.V()).value() << " "
        << "1/s" << endl;
}

void Foam::thermalHydraulicModels::onePhase::calcCumulContErr()
{
    if (pimple_.finalIter())
    {
        scalar& cumulContErr(fluid_.cumulContErr());

        scalarField integralContErr
        (
            mesh_.time().deltaTValue()*
            fvc::volumeIntegrate(fluid_.contErr())
        );

        const scalarField& V(mesh_.V());
        scalar totV(0);

        forAll(V, i)
        {
            cumulContErr += integralContErr[i];
            totV += V[i];
        }
        reduce(cumulContErr, sumOp<scalar>());
        reduce(totV, sumOp<scalar>());

        Info<< "Cumulative continuity error = " 
            << (cumulContErr/totV)
            << " kg/m3" << endl;
    }
}


// ************************************************************************* //