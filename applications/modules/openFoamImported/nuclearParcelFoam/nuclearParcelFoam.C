/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2512                                                  |
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

#include "nuclearParcelFoam.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(nuclearParcelFoam, 0);
    addToRunTimeSelectionTable
    (
        solver,
        nuclearParcelFoam,
        dynamicFvMesh
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::nuclearParcelFoam::nuclearParcelFoam
(
    dynamicFvMesh& mesh
)
:
    solver(mesh),
    pimple_(mesh_),
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
    Uf_
    (
        new surfaceVectorField
        (
            IOobject
            (
                "Uf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(U_)
        )
    ),
    laminarTransport_(U_, phi_),
    beta_
    (
        "beta",
        dimless/dimTemperature,
        laminarTransport_
    ),
    CpRef_
    (
        "CpRef",
        dimEnergy/dimMass/dimTemperature,
        laminarTransport_
    ),
    TRef_
    (
        "TRef",
        dimTemperature,
        laminarTransport_
    ),
    Pr_
    (
        "Pr",
        dimless,
        laminarTransport_
    ),
    Prt_
    (
        "Prt",
        dimless,
        laminarTransport_
    ),
    T_
    (
        IOobject
        (
            "T",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
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
    vort_
    (
        IOobject
        (
            "vort",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), vector::zero)
    ),
    turbulence_
    (
        incompressible::turbulenceModel::New(U_, phi_, laminarTransport_)
    ),
    rhoInfValue_
    (
        "rhoInf",
        dimDensity,
        laminarTransport_
    ),
    rhoInf_
    (
        IOobject
        (
            "rho",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        rhoInfValue_
    ),
    muc_
    (
        IOobject
        (
            "muc",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoInf_*laminarTransport_.nu()
    ),
    rhok_
    (
        IOobject
        (
            "rhok",
            mesh_.time().timeName(),
            mesh_
        ),
        1.0 - beta_*(T_ - TRef_)
    ),
    alphat_
    (
        IOobject
        (
            "alphat",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
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
        "ghf",
        (g_ & mesh_.Cf()) + mag(g_)*hRef_
    ),
    p_
    (
        IOobject
        (
            "p",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh_ + rhok_*gh_
    ),
    pRefCell_(0),
    pRefValue_(0.0),
    Qdot_
    (
        IOobject
        (
            "Qdot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
    ),
    MRF_(mesh_),
    radiation_(radiation::radiationModel::New(T_)),
    rhoCpRef_
    (
        "rhoCpRef",
        rhoInfValue_*CpRef_
    ),
    fvOptions_(fv::options::New(mesh_)),
    FP_
    (
        "FPCloud",
        rhoInf_,
        U_,
        muc_,
        vort_,
        T_,
        g_
    ),
    tsurfaceFilm_(regionModels::surfaceFilmModel::New(mesh_, g_, "surfaceFilm")),
    surfaceFilm_(tsurfaceFilm_()),
    cumulativeContErr_(0)
{
    turbulence_->validate();

    setRefCell
    (
        p_,
        p_rgh_,
        pimple_.dict(),
        pRefCell_,
        pRefValue_
    );

    if (p_rgh_.needReference())
    {
        p_ += dimensionedScalar
        (
            "p",
            p_.dimensions(),
            pRefValue_ - getRefCellValue(p_, pRefCell_)
        );
    }

    mesh_.setFluxRequired(p_rgh_.name());
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::nuclearParcelFoam::correctPhysics()
{
    FP_.storeGlobalPositions();
    FP_.evolve();

    surfaceFilm_.evolve();

    while (pimple_.loop())
    {
        #include "include/equations/UEqn.H"

        vort_ = fvc::curl(U_);
        #include "include/equations/TEqn.H"

        while (pimple_.correct())
        {
            #include "include/equations/pEqn.H"
        }

        if (pimple_.turbCorr())
        {
            laminarTransport_.correct();
            turbulence_->correct();
        }
    }
}


scalar Foam::solvers::nuclearParcelFoam::maxDeltaT()
{
    scalar newDeltaT = mesh_.time().controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

    scalar maxCo = mesh_.time().controlDict().getOrDefault<scalar>("maxCo", 1);

    if (mesh_.time().value() > mesh_.time().controlDict().get<scalar>("deltaT"))
    {
        scalar CoNum = 0.0;
        scalar meanCoNum = 0.0;

        if (mesh_.nInternalFaces())
        {
            scalarField sumPhi(fvc::surfaceSum(mag(phi_))().primitiveField());

            CoNum = 0.5*gMax(sumPhi/mesh_.V().field())*mesh_.time().deltaTValue();

            meanCoNum =
                0.5*(gSum(sumPhi)/gSum(mesh_.V().field()))*mesh_.time().deltaTValue();
        }

        Info<< "Courant Number mean: " << meanCoNum
            << " max: " << CoNum << endl;

        scalar maxDeltaTFact = maxCo/(CoNum + SMALL);
        scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

        newDeltaT = min(deltaTFact*mesh_.time().deltaTValue(), newDeltaT);
    }

    return newDeltaT;
}


// ************************************************************************* //
