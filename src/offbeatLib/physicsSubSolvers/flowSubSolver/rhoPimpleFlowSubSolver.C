/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "rhoPimpleFlowSubSolver.H"
#include "typeInfo.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fvc.H"
#include "CorrectPhi.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "hydrostaticInitialisation.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rhoPimpleFlowSubSolver, 0);
    addToRunTimeSelectionTable
    (
        flowSubSolver, 
        rhoPimpleFlowSubSolver, 
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


Foam::fvMesh& Foam::rhoPimpleFlowSubSolver::getMesh
(
    const Time& runTime
)
{
    word regionName = flowSubSolverDict_.lookup("region");
    
    if (runTime.foundObject<fvMesh>(regionName))
    {
        return runTime.lookupObjectRef<fvMesh>(regionName);
    }
    else
    {
        Info
            << "Create mesh for time = "
            << runTime.timeName() << nl << endl;
        
        meshPtr_.set(
            new fvMesh
            (
                IOobject
                (
                    regionName,
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );

        return meshPtr_();
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix>
Foam::rhoPimpleFlowSubSolver::solveMomentum()
{
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(rho_, U_) + fvm::div(phi_, U_)
      + turbulence_->divDevTau(U_)
    ==
        fvModels_.source(rho_, U_)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvConstraints_.constrain(UEqn);

    if (pimple_.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p_));

        fvConstraints_.constrain(U_);
        K_ = 0.5*magSqr(U_);
    }
    
    return tUEqn;
}


void Foam::rhoPimpleFlowSubSolver::correctPressure
(
    fvVectorMatrix& UEqn
)
{
    if ((!mesh_.steady() && !pimple_.simpleRho()) || pimple_.consistent())
    {
        rho_ = thermo_->rho();
    }

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution
    const volScalarField psip0(psi_*p_);

    const volScalarField rAU("rAU", 1.0/UEqn.A());
    const surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho_*rAU));

    tmp<volScalarField> rAtU
    (
        pimple_.consistent()
      ? volScalarField::New("rAtU", 1.0/(1.0/rAU - UEqn.H1()))
      : tmp<volScalarField>(nullptr)
    );
    tmp<surfaceScalarField> rhorAtUf
    (
        pimple_.consistent()
      ? surfaceScalarField::New("rhoRAtUf", fvc::interpolate(rho_*rAtU()))
      : tmp<surfaceScalarField>(nullptr)
    );

    const volScalarField& rAAtU = pimple_.consistent() ? rAtU() : rAU;
    const surfaceScalarField& rhorAAtUf =
        pimple_.consistent() ? rhorAtUf() : rhorAUf;

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U_, p_));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho_)*fvc::flux(HbyA)
      + rhorAUf*fvc::ddtCorr(rho_, U_, phi_)
    );

    fvc::makeRelative(phiHbyA, rho_, U_);

    bool adjustMass = false;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_, rho_, U_, phiHbyA, rhorAAtUf);

    if (mesh_.steady())
    {
        adjustMass = adjustPhi(phiHbyA, U_, p_);
    }

    if (pimple_.consistent())
    {
        phiHbyA += (rhorAAtUf - rhorAUf)*fvc::snGrad(p_)*mesh_.magSf();
        HbyA += (rAAtU - rAU)*fvc::grad(p_);
    }

    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(rho_) + psi_*correction(fvm::ddt(p_))
      + fvc::div(phiHbyA)
    ==
        fvModels_.source(psi_, p_, rho_.name())
    );

    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAAtUf, p_));

        pEqn.setReference
        (
            pressureReference_.refCell(),
            pressureReference_.refValue()
        );

        pEqn.solve();

        if (pimple_.finalNonOrthogonalIter())
        {
            phi_ = phiHbyA + pEqn.flux();
        }
    }

    continuityErrors();

    // Explicitly relax pressure for momentum corrector
    p_.relax();

    U_ = HbyA - rAAtU*fvc::grad(p_);
    U_.correctBoundaryConditions();
    fvConstraints_.constrain(U_);
    K_ = 0.5*magSqr(U_);

    fvConstraints_.constrain(p_);

    // For steady compressible closed-volume cases adjust the pressure level
    // to obey overall mass continuity
    if (adjustMass && !thermo_->incompressible())
    {
        p_ += (initialMass_ - fvc::domainIntegrate(thermo_->rho()))
            /fvc::domainIntegrate(psi_);
        p_.correctBoundaryConditions();
    }

    if (mesh_.steady() || pimple_.simpleRho() || adjustMass)
    {
        rho_ = thermo_->rho();
    }

    if (thermo_->dpdt())
    {
        dpdt_ = fvc::ddt(p_);
    }
}


void Foam::rhoPimpleFlowSubSolver::correctEnergy()
{
    //- Rename member variables for readibility
    volScalarField& he = thermo_->he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho_, he) + fvm::div(phi_, he)
      + fvc::ddt(rho_, K_) + fvc::div(phi_, K_)
      + (
            he.name() == "e"
          ? fvc::div(fvc::absolute(phi_, rho_, U_), p_/rho_)
          : -dpdt_
        )
      + thermoTransport_->divq(he)
     ==
        fvModels_.source(rho_, he)
    );

    EEqn.relax();

    fvConstraints_.constrain(EEqn);

    EEqn.solve();

    fvConstraints_.constrain(he);

    thermo_->correct();
}


void Foam::rhoPimpleFlowSubSolver::continuityErrors()
{
    volScalarField contErr(fvc::div(phi_));
 
    scalar sumLocalContErr = mesh_.time().deltaTValue()*
        mag(contErr)().weightedAverage(mesh_.V()).value();

    scalar globalContErr = mesh_.time().deltaTValue()*
        contErr.weightedAverage(mesh_.V()).value();
    cumulativeContErr_ += globalContErr;

    Info<< "time step continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr_
        << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rhoPimpleFlowSubSolver::rhoPimpleFlowSubSolver
(
    const Time& runTime,
    const dictionary& dict
)
:
    flowSubSolver(runTime, dict),
    meshPtr_(),
    mesh_(getMesh(runTime)),
    thermo_(fluidThermo::New(mesh_)),
    p_(thermo_->p()),
    psi_(thermo_->psi()),
    rho_(
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo_->rho()
    ),
    U_(
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    phi_(
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(rho_*U_) & mesh_.Sf()
    ),
    turbulence_(
        compressible::momentumTransportModel::New
        (
            rho_,
            U_,
            phi_,
            thermo_()
        )
    ),
    thermoTransport_
    (
        fluidThermophysicalTransportModel::New
        (
            turbulence_(),
            thermo_()
        )
    ),
    dpdt_
    (
        IOobject
        (
            "dpdt",
            runTime.timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(p_.dimensions()/dimTime, 0)
    ),
    K_("K", 0.5*magSqr(U_)),
    pimple_(mesh_),
    pressureReference_(p_, pimple_.dict()),
    initialMass_(fvc::domainIntegrate(rho_)),
    cumulativeContErr_(0),
    fvModels_(mesh_),
    fvConstraints_(mesh_)
{
    mesh_.setFluxRequired(p_.name());

    turbulence_->validate();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::rhoPimpleFlowSubSolver::init()
{
    Time& runTime = const_cast<Time&>(mesh_.time());
    
    pimple_.run(runTime);
}

void Foam::rhoPimpleFlowSubSolver::correct()
{
    {
        fvVectorMatrix UEqn = solveMomentum();
        
        while (pimple_.correct())
        {
            correctEnergy();
            correctPressure(UEqn);
        }
    }
    
    continuityErrors();

    if (pimple_.turbCorr())
    {
        turbulence_->correct();
        thermoTransport_->correct();
    }
}


void Foam::rhoPimpleFlowSubSolver::finalise()
{
    rho_ = thermo_->rho();

    if (!pimple_.turbCorr())
    {
        turbulence_->correct();
        thermoTransport_->correct();
    }
}


bool Foam::rhoPimpleFlowSubSolver::converged()
{
    return not pimple_.loop();
}


Foam::scalar Foam::rhoPimpleFlowSubSolver::nextDeltaT() const
{
    if (pimple_.dict().found("maxCo"))
    {
        scalar maxCo = pimple_.dict().lookup<scalar>("maxCo");
        
        notImplemented("rhoPimpleFlowSubSolver::nextDeltaT() for specified maxCo");
        
        return great;
    }
    else
    {
        return great;
    }
}

// ************************************************************************* //
