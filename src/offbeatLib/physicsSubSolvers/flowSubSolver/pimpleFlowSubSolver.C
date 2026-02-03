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

#include "pimpleFlowSubSolver.H"
#include "typeInfo.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fvc.H"
#include "CorrectPhi.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pimpleFlowSubSolver, 0);
    addToRunTimeSelectionTable
    (
        flowSubSolver, 
        pimpleFlowSubSolver, 
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


Foam::fvMesh& Foam::pimpleFlowSubSolver::getMesh
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
Foam::pimpleFlowSubSolver::solveMomentum()
{
    volVectorField& U = U_;
    volScalarField& p = p_;
    surfaceScalarField& phi = phi_;
    
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + turbulence_->divDevSigma(U)
     ==
        fvModels().source(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvConstraints().constrain(UEqn);

    if (pimple_.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvConstraints().constrain(U);
    }
    
    return tUEqn;
}


void Foam::pimpleFlowSubSolver::correctPressure
(
    fvVectorMatrix& UEqn
)
{
    volVectorField& U = U_;
    volScalarField& p = p_;
    surfaceScalarField& phi = phi_;
    
    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
    );

    if (p.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p);
        fvc::makeAbsolute(phiHbyA, U);
    }
    
    tmp<volScalarField> rAtU(rAU);

    if (pimple_.consistent())
    {
        rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
        phiHbyA +=
            fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh_.magSf();
        HbyA -= (rAU - rAtU())*fvc::grad(p);
    }

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAtU());

    // Non-orthogonal pressure corrector loop
    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
        );

        pEqn.setReference
        (
            pressureReference_.refCell(),
            pressureReference_.refValue()
        );

        pEqn.solve();

        if (pimple_.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    // Explicitly relax pressure for momentum corrector
    p.relax();

    U = HbyA - rAtU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvConstraints().constrain(U);

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);
}


void Foam::pimpleFlowSubSolver::continuityErrors()
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

Foam::pimpleFlowSubSolver::pimpleFlowSubSolver
(
    const Time& runTime,
    const dictionary& dict
)
:
    flowSubSolver(runTime, dict),
    meshPtr_(),
    mesh_(getMesh(runTime)),
    p_(
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
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
        fvc::flux(U_)
    ),
    laminarTransport_(U_, phi_),
    turbulence_(
        incompressible::momentumTransportModel::New
        (
            U_, 
            phi_, 
            laminarTransport_
        )
    ),
    pimple_(mesh_),
    pressureReference_(p_, pimple_.dict()),
    cumulativeContErr_(0),
    fvModels_(mesh_),
    fvConstraints_(mesh_)
{
    mesh_.setFluxRequired(p_.name());

    turbulence_->validate();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::pimpleFlowSubSolver::init()
{
    Time& runTime = const_cast<Time&>(mesh_.time());
    
    pimple_.run(runTime);
}

void Foam::pimpleFlowSubSolver::correct()
{
    {
        fvVectorMatrix UEqn = solveMomentum();
        
        while (pimple_.correct())
        {
            correctPressure(UEqn);
        }
    }
    
    continuityErrors();

    if (pimple_.turbCorr())
    {
        laminarTransport_.correct();
        turbulence_->correct();
    }
}


void Foam::pimpleFlowSubSolver::finalise()
{
    if (not pimple_.turbCorr())
    {
        laminarTransport_.correct();
        turbulence_->correct();
    }
}


bool Foam::pimpleFlowSubSolver::converged()
{
    return not pimple_.loop();
}


Foam::scalar Foam::pimpleFlowSubSolver::nextDeltaT() const
{
    // TODO
    return great;
}

// ************************************************************************* //
