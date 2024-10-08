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

#include "solver.H"
#include "localEulerDdtScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solver, 0);
    defineRunTimeSelectionTable(solver, fvMesh);
}


Foam::scalar Foam::solver::deltaTFactor = 1.2;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solver::writeData(Ostream&) const
{
    NotImplemented;
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solver::solver(fvMesh& mesh)
:
    regIOobject
    (
        IOobject
        (
            typeName,
            mesh.time().timeName(),
            mesh
        )
    ),
    mesh_(mesh),
    steady(mesh_.schemes().steady()),
    LTS(fv::localEulerDdt::enabled(mesh)),
    // fvModelsPtr(nullptr),
    // fvConstraintsPtr(nullptr),
    mesh(mesh_),
    runTime(mesh_.time())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solver::~solver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Foam::fvModels& Foam::solver::fvModels() const
// {
//     if (!fvModelsPtr)
//     {
//         fvModelsPtr = &Foam::fvModels::New(mesh);
//     }

//     return *fvModelsPtr;
// }


// Foam::fvConstraints& Foam::solver::fvConstraints() const
// {
//     if (!fvConstraintsPtr)
//     {
//         fvConstraintsPtr = &Foam::fvConstraints::New(mesh);
//     }

//     return *fvConstraintsPtr;
// }


// ************************************************************************* //
