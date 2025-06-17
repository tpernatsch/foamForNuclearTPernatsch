/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2406                                                  |
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

#include "compressibleInterFoam.H"
#include "addToRunTimeSelectionTable.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(compressibleInterFoam, 0);
    addToRunTimeSelectionTable
    (
        solver, 
        compressibleInterFoam, 
        dynamicFvMesh
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::compressibleInterFoam::compressibleInterFoam
(
    dynamicFvMesh& mesh_
)
:
    solver(mesh_),
    mesh_(mesh_),
    pimple_(mesh_),
    p_rgh_
    (
        IOobject
        (
            "p_rgh",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh_
    ),
    U_
    (
        IOobject
        (
            "U",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    phi_
    (
        IOobject
        (
            "phi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U_)
    ),
    mixture_(U_,phi_),
    alpha1_(mixture_.alpha1()),
    alpha2_(mixture_.alpha2()),
    rho1_(mixture_.thermo1().rho()),
    rho2_(mixture_.thermo2().rho()),
    rho_
    (
        IOobject
        (
            "rho",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha1_*rho1_ + alpha2_*rho2_
    ),
    pMin_
    (
        "pMin",
        dimPressure,
        mixture_
    ),
    g_
    (
        IOobject
        (
            "g",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    hRef_
    (
        IOobject
        (
            "hRef",
            mesh_.time().constant(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar("hRef", dimLength, 0)
    ),
    gh_
    (
        "gh",
        (g_ & mesh_.C()) + mag(g_)*hRef_
    ),
    ghf_
    (
        "ghf_",
        (g_ & mesh_.Cf()) + mag(g_)*hRef_
    ),
    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho_)*phi_
    ),
    dgdt_
    (
        "dgdt",
        alpha1_*fvc::div(phi_)
    ),
    K_
    (
        "K",
        0.5*magSqr(U_)
    ),
    fvOptions_(fv::options::New(mesh_)),
    p_(mixture_.p()),
    T_(mixture_.T()),
    psi1_(mixture_.thermo1().psi()),
    psi2_(mixture_.thermo2().psi()),
    alphaPhi10Header
    (
        IOobject::groupName("alphaPhi0", alpha1_.group()),
        mesh_.time().timeName(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alphaRestart( alphaPhi10Header.typeHeaderOk<surfaceScalarField>(true)),
    alphaPhi10
    (
        alphaPhi10Header,
        phi_*fvc::interpolate(alpha1_)
    ),
    alphaPhiUn
    (
        IOobject
        (
            "alphaPhiUn",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(phi_.dimensions(), Zero)
    ),
    MRF_(mesh_)
{

    turbulence_.reset(new compressibleInterPhaseTransportModel
    (
        rho_,
        U_,
        phi_,
        rhoPhi_,
        alphaPhi10,
        mixture_ 
    ));

    Info << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- Solve according to flags
void Foam::solvers::compressibleInterFoam::correctPhysics()
{   
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple_.loop())
        {
            #include "alphaControls.H"
            #include "compressibleAlphaEqnSubCycle.H"

            turbulence_->correctPhasePhi();

            #include "UEqn.H"
            volScalarField divUp("divUp", fvc::div(fvc::absolute(phi_, U_), p_));
            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple_.correct())
            {
                #include "pEqn.H"
            }

            if (pimple_.turbCorr())
            {
                turbulence_->correct();
            }
        }
}

void Foam::solvers::compressibleInterFoam::correctTightlyCoupledPhysics()
{

}

void Foam::solvers::compressibleInterFoam::correctFluidMechanics()
{
}

void Foam::solvers::compressibleInterFoam::correctEnergy()
{
}



void Foam::solvers::compressibleInterFoam::correctCourant()
{

}

void Foam::solvers::compressibleInterFoam::correctContErr()
{

}


void Foam::solvers::compressibleInterFoam::printContErr()
{

}

void Foam::solvers::compressibleInterFoam::calcCumulContErr()
{

}


scalar Foam::solvers::compressibleInterFoam::maxDeltaT()
{
    return scalar(VGREAT);
}


// ************************************************************************* //