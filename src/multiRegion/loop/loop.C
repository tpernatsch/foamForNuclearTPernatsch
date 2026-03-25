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

#include "loop.H"
#include "zeroGradientFvPatchFields.H"
#include "mergeOrSplitBaffles.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"
#include "fixedValuePointPatchFields.H"
#include "primitivePatchInterpolation.H"
#include "mappedPatchBase.H"
#include "tractionDisplacementFvPatchVectorField.H"
#include "typeInfo.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(loop, 0);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::loop::loop
(
    dynamicFvMesh& mesh
)
:
    solver(mesh),
    IOdictionary
    (
        IOobject
        (
            "regionsDict",
            runTime().time().system(),
            runTime().db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    meshHandler_(nullptr)
{
}


// * * * * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * * * //

void Foam::solvers::loop::correctBaffleLessFields()
{
    forAll(solvers_, solvI)
    {
        solvers_[solvI].correctBaffleLessFields();
    }
}


void Foam::solvers::loop::deformMesh()
{
    forAll(solvers_, solvI)
    {
        solvers_[solvI].deformMesh();
    }
}


scalar Foam::solvers::loop::getResidual()
{
    scalar residual = 0;
    forAll(solvers_, solvI)
    {
        residual = max(residual, solvers_[solvI].getResidual());
    }

    return residual;
}

scalar Foam::solvers::loop::maxDeltaT()
{
    scalar maxDeltaT= VGREAT;
    forAll(solvers_, solvI)
    {
        maxDeltaT = min(maxDeltaT, solvers_[solvI].maxDeltaT());
    }

    return maxDeltaT;
}


// ************************************************************************* //
