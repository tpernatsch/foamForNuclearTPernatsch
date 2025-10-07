/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2506                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2025 OpenCFD Ltd.          |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "onePhase.H"
#include "addToRunTimeSelectionTable.H"
#include "regimeMapModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(onePhase, 0);
    addToRunTimeSelectionTable
    (
        solver,
        onePhase,
        dynamicFvMesh
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::onePhase::onePhase
(
    dynamicFvMesh& mesh
)
:
    thermalHydraulicsModel
    (
        mesh.time(),
        mesh,
        fv::options::New(mesh)
    ),
    //- The structure needs to be created before the fluids because
    //  the turbulence models created by the fluids might require a reference
    //  to the structure (obtained via objectRegistry lookup in the specific
    //  turbulence model class)
    structure_
    (
        this->subDict("structureProperties"),
        mesh,
        this->powerDensityNeutronics_
    ),
    fluid_
    (
        (this->found("fluidProperties"))
    ?   this->subDict("fluidProperties") : *this,
        mesh,
        word(""),   //- This is the phase name, setting it to "" signals a
                    //  onePhase solver to the rest of the FFS library
        this->powerDensityNeutronicsToLiquid_,
        false       //- No need to read or write the fluid phaseFraction in
                    //  onePhase, it is tied to the structure phaseFraction
    ),
    FSPair_(fluid_, structure_, *this),
    residual_(0)
{

    //- Create turbulence model
    fluid_.constructTurbulenceModel();

    //- Set phase fraction fields (constant in time), structure has priority
    fluid_.volScalarField::operator=(1.0-structure_);

    //- The normalized field is non-trivial (i.e. different than 1) only in the
    //  twoPhase solver. However, it is used by some models in the shared
    //  thermal-hydraulics library, so it should be set nonetheless! The most
    //  important quantity that relies on this is the Reynolds computed by
    //  the fluidStructurePair object
    fluid_.normalized() = fluid_/(1.0-structure_);

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

    //- Compute initialFluidMass
    initialFluidMass_ = fvc::domainIntegrate(fluid_.rho()*fluid_);

    Info << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- Solve according to flags
void Foam::solvers::onePhase::correctPhysics()
{

    residual_=0;

    bool solveEnergy(mesh_.solutionDict().subDict("PIMPLE").get<bool>("solveEnergy"));
    bool solveFluidMechanics(mesh_.solutionDict().subDict("PIMPLE").get<bool>("solveFluidMechanics"));

    Info << "Region :" << mesh_.name()<<nl<<endl;


    while(pimple_.loop())
    {

        if (solveFluidMechanics or solveEnergy)
            correctModels(true, true);

        if(solveFluidMechanics)
            correctFluidMechanics();

        if(solveEnergy)
            correctEnergy();

        Info << endl;
    }
}

void Foam::solvers::onePhase::correctTightlyCoupledPhysics()
{
    Info <<nl;
    bool solveEnergy(mesh_.solutionDict().subDict("PIMPLE").get<bool>("solveEnergy"));
    correctModels(true,true);

    if(solveEnergy)
        correctEnergy();

    Info <<nl;
}

void Foam::solvers::onePhase::correctFluidMechanics()
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

void Foam::solvers::onePhase::correctEnergy()
{

    #include "EEqn_1p.H"
}

void Foam::solvers::onePhase::correctModels
(
    bool solveFluidDynamics,
    bool solveEnergy
)
{
    this->correctRegimeMaps();
    FSPair_.correct(solveFluidDynamics, solveEnergy);
}

void Foam::solvers::onePhase::correctCourant()
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

void Foam::solvers::onePhase::correctContErr()
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


void Foam::solvers::onePhase::printContErr()
{
    volScalarField contErrRel(fluid_.contErr()/fluid_.rho());

    Info<< "Instantaneous relative continuity error (avg) = "
        << contErrRel.weightedAverage(mesh_.V()).value()
        //<< " " << min(contErrRel).value()
        //<< " " << max(contErrRel).value()
        << " 1/s" << endl;
}

void Foam::solvers::onePhase::calcCumulContErr()
{
    if (pimple_.finalIter())
    {
        scalar& cumulContErr(fluid_.cumulContErr());
        const volScalarField& cE(fluid_.contErr());

        const scalarField& V(mesh_.V());
        const scalar& dt(mesh_.time().deltaTValue());

        scalar totV(0);
        scalar deltaCumulContErr(0);

        forAll(V, i)
        {
            const scalar& Vi(V[i]);
            deltaCumulContErr += cE[i]*Vi*dt;
            totV += Vi;
        }
        reduce(deltaCumulContErr, sumOp<scalar>());
        reduce(totV, sumOp<scalar>());

        cumulContErr += deltaCumulContErr;

        Info<< "Cumulative continuity error = "
            << (cumulContErr/totV)
            << " kg/m3" << endl;
    }
}


scalar Foam::solvers::onePhase::maxDeltaT()
{

    this->correctCourant();
    scalar newDeltaT = runTime_.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);
    scalar maxCo =runTime_.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
    scalar maxDeltaTFact = maxCo/(CoNum_ + SMALL);
    scalar deltaTFact =  min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);
    return min(deltaTFact*runTime_.deltaTValue(),newDeltaT);
}


// ************************************************************************* //