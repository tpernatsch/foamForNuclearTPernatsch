/*---------------------------------------------------------------------------*\
|                                                              |
|      / /      / | / /          / /               |
|     / /    /  \  /  |/ /    / /     /  \ /  `/  /  ` \    |
|    / // /  /  / / /|  /  // / /    / // // // /  / / / / / /    |
|    \/   \/ // |/          //       \/ \,/  // // //     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2412                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2024 OpenCFD Ltd.          |
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

#include "compressibleInterDyMFoam.H"
#include "addToRunTimeSelectionTable.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(compressibleInterDyMFoam, 0);
    addToRunTimeSelectionTable
    (
        solver,
        compressibleInterDyMFoam,
        dynamicFvMesh
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::compressibleInterDyMFoam::compressibleInterDyMFoam
(
    dynamicFvMesh& mesh
)
:
    solver(mesh),
    dynMesh(mesh),
    pimple(dynMesh),
    cumulativeContErr(0),
    p_rgh
    (
        IOobject
        (
            "p_rgh",
            dynMesh.time().timeName(),
            dynMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        dynMesh
    ),
    U
    (
        IOobject
        (
            "U",
            dynMesh.time().timeName(),
            dynMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        dynMesh
    ),
    phi
    (
        IOobject
        (
            "phi",
            dynMesh.time().timeName(),
            dynMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U)
    ),
    mixture(U,phi),
    alpha1(mixture.alpha1()),
    alpha2(mixture.alpha2()),
    rho1(mixture.thermo1().rho()),
    rho2(mixture.thermo2().rho()),
    rho
    (
        IOobject
        (
            "rho",
            dynMesh.time().timeName(),
            dynMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha1*rho1 + alpha2*rho2
    ),
    pMin
    (
        "pMin",
        dimPressure,
        mixture
    ),
    g
    (
        IOobject
        (
            "g",
            dynMesh.time().constant(),
            dynMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    hRef
    (
        IOobject
        (
            "hRef",
            dynMesh.time().constant(),
            dynMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar("hRef", dimLength, 0)
    ),
    ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRef
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    ),
    gh
    (
        "gh",
        (g & dynMesh.C()) + mag(g)*hRef
    ),
    ghf
    (
        "ghf",
        (g & dynMesh.Cf()) + mag(g)*hRef
    ),
    rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            dynMesh.time().timeName(),
            dynMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)*phi
    ),
    dgdt
    (
        "dgdt",
        alpha1*fvc::div(phi)
    ),
    Uf
    (
        IOobject
        (
            "Uf",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U)
    ),
    K
    (
        "K",
        0.5*magSqr(U)
    ),
    fvOptions(fv::options::New(dynMesh)),
    p(mixture.p()),
    T(mixture.T()),
    psi1(mixture.thermo1().psi()),
    psi2(mixture.thermo2().psi()),
    alphaPhi10Header
    (
        IOobject::groupName("alphaPhi0", alpha1.group()),
        dynMesh.time().timeName(),
        dynMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alphaRestart( alphaPhi10Header.typeHeaderOk<surfaceScalarField>(true)),
    alphaPhi10
    (
        alphaPhi10Header,
        phi*fvc::interpolate(alpha1)
    ),
    alphaPhiUn
    (
        IOobject
        (
            "alphaPhiUn",
            dynMesh.time().timeName(),
            dynMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dynMesh,
        dimensionedScalar(phi.dimensions(), Zero)
    ),
    MRF(dynMesh)
{

    turbulence.reset(new compressibleInterPhaseTransportModel
    (
        rho,
        U,
        phi,
        rhoPhi,
        alphaPhi10,
        mixture
    ));

    if (U.nOldTimes())
    {
        volVectorField* Uold = &U.oldTime();
        volScalarField* Kold = &K.oldTime();
        *Kold == 0.5*magSqr(*Uold);
    
        while (Uold->nOldTimes())
        {
            Uold = &Uold->oldTime();
            Kold = &Kold->oldTime();
            *Kold == 0.5*magSqr(*Uold);
        }
    }


    Info << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- Solve according to flags
void Foam::solvers::compressibleInterDyMFoam::correctPhysics()
{


        // --- Pressure-velocity pimple corrector loop

    bool correctPhi
    (
        pimple.dict().getOrDefault("correctPhi", true)
    );
    
    bool checkMeshCourantNo
    (
        pimple.dict().getOrDefault("checkMeshCourantNo", false)
    );
    
    bool moveMeshOuterCorrectors
    (
        pimple.dict().getOrDefault("moveMeshOuterCorrectors", false)
    );

    // Store divU and divUp from the previous dynMesh so that it can be mapped
    // and used in correctPhi to ensure the corrected phi has the
    // same divergence
    volScalarField divU("divU0", fvc::div(fvc::absolute(phi, U)));
    volScalarField divUp("divUp", fvc::div(fvc::absolute(phi, U), p));

    // --- Pressure-velocity pimple corrector loop
    while (pimple.loop())
    {
        if (pimple.firstIter() || moveMeshOuterCorrectors)
        {
            scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

            dynMesh.update();

            if (dynMesh.changing())
            {
                MRF.update();

                Info<< "Execution time for dynMesh.update() = "
                    << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                    << " s" << endl;

                gh = (g & dynMesh.C()) - ghRef;
                ghf = (g & dynMesh.Cf()) - ghRef;
            }

            if ((dynMesh.changing() && correctPhi))
            {
                // Calculate absolute flux from the mapped surface velocity
                phi = dynMesh.Sf() & Uf;

                #include "correctPhi.H"

                // Make the fluxes relative to the dynMesh motion
                fvc::makeRelative(phi, U);

                mixture.correct();
            }

            if (dynMesh.changing() && checkMeshCourantNo)
            {
                #include "meshCourantNo.H"
            }
        }

        #include "alphaControls.H"
        #include "compressibleAlphaEqnSubCycle.H"

        turbulence().correctPhasePhi();

        #include "UEqn.H"
        #include "TEqn.H"

        // --- Pressure corrector loop
        while (pimple.correct())
        {
            #include "pEqn.H"
        }

        if (pimple.turbCorr())
        {
            turbulence().correct();
        }
    }

    rho = alpha1*rho1 + alpha2*rho2;

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh = p - rho*gh;
    p_rgh.correctBoundaryConditions();
}

void Foam::solvers::compressibleInterDyMFoam::correctTightlyCoupledPhysics()
{
}

void Foam::solvers::compressibleInterDyMFoam::correctFluidMechanics()
{
}

void Foam::solvers::compressibleInterDyMFoam::correctEnergy()
{
}

void Foam::solvers::compressibleInterDyMFoam::correctCourant()
{
}

void Foam::solvers::compressibleInterDyMFoam::correctContErr()
{
}

void Foam::solvers::compressibleInterDyMFoam::printContErr()
{
}

void Foam::solvers::compressibleInterDyMFoam::calcCumulContErr()
{
}

scalar Foam::solvers::compressibleInterDyMFoam::maxDeltaT()
{
    volScalarField divU("divU0", fvc::div(fvc::absolute(phi, U)));
    volScalarField divUp("divUp", fvc::div(fvc::absolute(phi, U), p));

    bool adjustTimeStep(false);
    scalar maxCo(1);
    scalar maxDeltaT(GREAT);

    bool correctPhi
    (
        pimple.dict().getOrDefault("correctPhi", true)
    );
    
    bool checkMeshCourantNo
    (
        pimple.dict().getOrDefault("checkMeshCourantNo", false)
    );
    
    bool moveMeshOuterCorrectors
    (
        pimple.dict().getOrDefault("moveMeshOuterCorrectors", false)
    );

    #include "readDyMControls.H"

    if (LTS)
    {
        Info <<"Ciao fratm" << endl;
    }

    #include "CourantNo.H"
    #include "alphaCourantNo.H"
    scalar maxDeltaTFact =
    min(maxCo/(CoNum + SMALL), maxAlphaCo/(alphaCoNum + SMALL));

    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);
        
    return (min
        (
            deltaTFact*runTime.deltaTValue(),
            maxDeltaT
        ));
}


// ************************************************************************* //