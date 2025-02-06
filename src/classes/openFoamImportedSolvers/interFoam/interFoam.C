/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
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

#include "interFoam.H"
#include "addToRunTimeSelectionTable.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(interFoam, 0);
    addToRunTimeSelectionTable
    (
        solver,
        interFoam,
        fvMesh
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::interFoam::interFoam
(
    fvMesh& mesh
)
:
    solver(mesh),
    mesh(static_cast<dynamicFvMesh&>(mesh)),
    pimple(mesh),
    cumulativeContErr(0),
    p_rgh
    (
        IOobject
        (
            "p_rgh",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    U
    (
        IOobject
        (
            "U",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    phi
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U)
    ),
    mixture(U,phi),
    alpha1(mixture.alpha1()),
    alpha2(mixture.alpha2()),
    rho1(mixture.rho1()),
    rho2(mixture.rho2()),
    rho
    (
        IOobject
        (
            "rho",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha1*rho1 + alpha2*rho2
    ),
    rhoPhi
    (
        IOobject
        (
            "rho",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(rho)*phi
    ),
    g
    (
        IOobject
        (
            "g",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    hRef
    (
        IOobject
        (
            "hRef",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(dimLength, Zero)  
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
        (g & mesh.C()) + mag(g)*hRef
    ),
    ghf
    (
        "ghf",
        (g & mesh.Cf()) + mag(g)*hRef
    ),
    p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*gh
    ),
    pRefCell(0),
    pRefValue(0),
    alphaPhi10Header
    (
        IOobject::groupName("alphaPhi0", alpha1.group()),
        mesh.time().timeName(),
        mesh,
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
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(phi.dimensions(), Zero)
    ),
    MRF(mesh),
    fvOptions(fv::options::New(mesh))
{


    bool correctPhi
    (
        pimple.dict().getOrDefault("correctPhi", true)
    );
    

    if (correctPhi)
    {
        rAU = new volScalarField
        (
            IOobject
            (
                "rAU",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("rAU", dimTime/dimDensity, 1)
        );

        #include "correctPhi.H"
    }
    else
    {
        CorrectPhi
        (
            U,
            phi,
            p_rgh,
            dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1),
            geometricZeroField(),
            pimple
        );

        {
            volScalarField contErr(fvc::div(phi));
        
            scalar sumLocalContErr = runTime.deltaTValue()*
                mag(contErr)().weightedAverage(mesh.V()).value();
        
            scalar globalContErr = runTime.deltaTValue()*
                contErr.weightedAverage(mesh.V()).value();
            cumulativeContErr += globalContErr;
        
            Info<< "time step continuity errors : sum local = " << sumLocalContErr
                << ", global = " << globalContErr
                << ", cumulative = " << cumulativeContErr
                << endl;
        }
    }


    turbulence.reset(new transportModelType
    (
        rho,
        U,
        phi,
        rhoPhi,
        mixture
    ));

    Info << endl;

    if (mesh.dynamic())
    {
        Info<< "Constructing face velocity Uf\n" << endl;
    
        Uf.reset
        (
            new surfaceVectorField
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
            )
        );
    }

    setRefCell
    (
        p,
        p_rgh,
        pimple.dict(),
        pRefCell,
        pRefValue
    );

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- Solve according to flags
void Foam::solvers::interFoam::correctPhysics()
{
        // --- Pressure-velocity PIMPLE corrector loop

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
    
    
    bool massFluxInterpolation
    (
        pimple.dict().getOrDefault("massFluxInterpolation", false)
    );
    
    bool adjustFringe
    (
        pimple.dict().getOrDefault("oversetAdjustPhi", false)
    );
    
    bool ddtCorr
    (
        pimple.dict().getOrDefault("ddtCorr", true)
    );
    while (pimple.loop())
    {
        if (pimple.firstIter() || moveMeshOuterCorrectors)
        {
            mesh.update();

            if (mesh.changing())
            {
                // Do not apply previous time-step mesh compression flux
                // if the mesh topology changed
                if (mesh.topoChanging())
                {
                    talphaPhi1Corr0.clear();
                }

                gh = (g & mesh.C()) - ghRef;
                ghf = (g & mesh.Cf()) - ghRef;

                MRF.update();

                if (correctPhi)
                {
                    // Calculate absolute flux
                    // from the mapped surface velocity
                    phi = mesh.Sf() & Uf();

                    #include "correctPhi.H"

                    // Make the flux relative to the mesh motion
                    fvc::makeRelative(phi, U);

                    mixture.correct();
                }

                if (checkMeshCourantNo)
                {
                    #include "meshCourantNo.H"
                }
            }
        }

        #include "alphaControls.H"
        #include "alphaEqnSubCycle.H"

        mixture.correct();

        if (pimple.frozenFlow())
        {
            continue;
        }

        #include "UEqn.H"

        // --- Pressure corrector loop
        while (pimple.correct())
        {
            #include "pEqn.H"
        }

        if (pimple.turbCorr())
        {
            turbulence->correct();
        }
    }
    // if(mesh.time().outputTime())
    //     mesh.write();
}

void Foam::solvers::interFoam::correctTightlyCoupledPhysics()
{
}

void Foam::solvers::interFoam::correctFluidMechanics()
{
}

void Foam::solvers::interFoam::correctEnergy()
{
}

void Foam::solvers::interFoam::correctCourant()
{
}

void Foam::solvers::interFoam::correctContErr()
{
}

void Foam::solvers::interFoam::printContErr()
{
}

void Foam::solvers::interFoam::calcCumulContErr()
{
}

scalar Foam::solvers::interFoam::maxDeltaT()
{
    scalar newDeltaT = mesh.time().controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

    // adjustTimeStep = runTime.controlDict().getOrDefault("adjustTimeStep", false);

    scalar maxCo = mesh.time().controlDict().getOrDefault<scalar>("maxCo", 1);

    if (mesh.time().value() > mesh.time().controlDict().get<scalar>("deltaT"))
    {
        scalar CoNum = 0.0;
        scalar meanCoNum = 0.0;

        if (mesh.nInternalFaces())
        {
            scalarField sumPhi
            (
                mixture.nearInterface()().primitiveField()
                *fvc::surfaceSum(mag(phi))().primitiveField()
            );

            CoNum = 0.5*gMax(sumPhi/mesh.V().field())*mesh.time().deltaTValue();

            meanCoNum =
                0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*mesh.time().deltaTValue();
        }

        Info<< "Courant Number mean: " << meanCoNum
            << " max: " << CoNum << endl;


        scalar maxAlphaCo
        (
            mesh.time().controlDict().get<scalar>("maxAlphaCo")
        );

        scalar alphaCoNum = 0.0;
        scalar meanAlphaCoNum = 0.0;

        if (mesh.nInternalFaces())
        {
            scalarField sumPhi
            (
                mixture.nearInterface()().primitiveField()
                *fvc::surfaceSum(mag(phi))().primitiveField()
            );

            alphaCoNum = 0.5*gMax(sumPhi/mesh.V().field())*mesh.time().deltaTValue();

            meanAlphaCoNum =
                0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*mesh.time().deltaTValue();
        }

        Info<< "Interface Courant Number mean: " << meanAlphaCoNum
            << " max: " << alphaCoNum << endl;

        scalar maxDeltaTFact =
            min(maxCo/(CoNum + SMALL), maxAlphaCo/(alphaCoNum + SMALL));

        scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

        newDeltaT =
        (
            min
            (
                deltaTFact*mesh.time().deltaTValue(),
                newDeltaT
            )
        );
    }

    return newDeltaT;
}


// ************************************************************************* //