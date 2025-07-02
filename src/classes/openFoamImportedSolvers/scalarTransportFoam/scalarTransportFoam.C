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

#include "scalarTransportFoam.H"
#include "addToRunTimeSelectionTable.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(scalarTransportFoam, 0);
    addToRunTimeSelectionTable
    (
        solver,
        scalarTransportFoam,
        fvMesh
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::scalarTransportFoam::scalarTransportFoam
(
    fvMesh& mesh_
)
:
    solver(mesh_),
    mesh_(mesh_),
    pimple_(mesh_),
    T_
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
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
    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    DT_("DT", dimViscosity, transportProperties_),
    fvOptions_(fv::options::New(mesh_))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::scalarTransportFoam::correctPhysics()
{
    phi_ = fvc::flux(U_);

    while (pimple_.loop())
    {
        while (pimple_.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T_)
              + fvm::div(phi_, T_)
              - fvm::laplacian(DT_, T_)
             ==
                fvOptions_(T_)
            );

            TEqn.relax();
            fvOptions_.constrain(TEqn);
            TEqn.solve();
            fvOptions_.correct(T_);
        }
    }
}

void Foam::solvers::scalarTransportFoam::correctTightlyCoupledPhysics()
{
}

void Foam::solvers::scalarTransportFoam::correctFluidMechanics()
{
}

void Foam::solvers::scalarTransportFoam::correctEnergy()
{
}

void Foam::solvers::scalarTransportFoam::correctCourant()
{
}

void Foam::solvers::scalarTransportFoam::correctContErr()
{
}

void Foam::solvers::scalarTransportFoam::printContErr()
{
}

void Foam::solvers::scalarTransportFoam::calcCumulContErr()
{
}

scalar Foam::solvers::scalarTransportFoam::maxDeltaT()
{
    scalar newDeltaT = mesh_.time().controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

    return newDeltaT;
}


// ************************************************************************* //