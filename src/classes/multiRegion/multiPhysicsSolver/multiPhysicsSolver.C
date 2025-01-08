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

#include "multiPhysicsSolver.H"
#include "zeroGradientFvPatchFields.H"
#include "mergeOrSplitBaffles.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(multiPhysicsSolver, 0);
    addToRunTimeSelectionTable
    (
        solver,
        multiPhysicsSolver,
        fvMesh
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::multiPhysicsSolver::multiPhysicsSolver
(
    fvMesh& mesh
)
:
    solver(mesh),
    IOdictionary
    (
        IOobject
        (
            "multiRegionCouplingDict",
            runTime.time().constant(),
            runTime.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    meshHandler_(nullptr)
{
}


// * * * * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * * * //

void Foam::solvers::multiPhysicsSolver::createSolvers(word name)
{
    multiPhysicsDict_ = this->subDict("multiPhysicsSolvers").subDict(name);
    solverNames_ = multiPhysicsDict_.subDict("solvers").toc();
    singlePhysicsSolverNames_.setSize(0);
    minResidual_ = multiPhysicsDict_.get<scalar>("minResidual");
    maxIterations_ = multiPhysicsDict_.get<label>("maxIterations");

    // First create the new solvers

    solvers_.setSize(solverNames_.size());

    correctOnlyEnergy_.setSize(solverNames_.size());
    forAll(correctOnlyEnergy_, boolI)
    {
        correctOnlyEnergy_[boolI]=true;
    }

    wordList energyOnlyList = multiPhysicsDict_.getOrDefault<wordList>("includeFluidMechanicsInLoop", wordList());

    forAll(solverNames_, nameI)
    {

        if(meshHandler_->contains(energyOnlyList, solverNames_[nameI]))
            correctOnlyEnergy_[nameI]=false;

        Info<< "Creating sub-scale solver " << solverNames_[nameI] << nl
            << endl;

        word solverType
        (
            multiPhysicsDict_
                .subDict("solvers")
                .get<word>(solverNames_[nameI])
        );

        // This replicates the structure of regionSolvers. This allows to create
        // subSolvers which are multiPhysicsSolvers themselves and so on,
        // creating a tree-like structure

        if (solverType != "multiPhysicsSolver")
        {
            solvers_.set
            (
                nameI,
                solver::New(solverType, meshHandler_->returnMesh(solverNames_[nameI]))
            );
            singlePhysicsSolverNames_.append(solverNames_[nameI]);
        }
        else
        {
            solvers_.set
            (
                nameI,
                solver::New(solverType, meshHandler_->returnMesh("dummy"))
            );
            Foam::solvers::multiPhysicsSolver* multiPhysicsSolverPtr = dynamic_cast<Foam::solvers::multiPhysicsSolver*>(&solvers_[nameI]);
            multiPhysicsSolverPtr->getMapper(*meshHandler_);
            multiPhysicsSolverPtr->createSolvers(solverNames_[nameI]);
            multiPhysicsSolverPtr = nullptr;
        }
    }
}


void Foam::solvers::multiPhysicsSolver::correctBaffleLessFields()
{
    forAll(solvers_, solvI)
    {
        solvers_[solvI].correctBaffleLessFields();
    }
}


void Foam::solvers::multiPhysicsSolver::deformMesh()
{
    forAll(solvers_, solvI)
    {
        solvers_[solvI].deformMesh();
    }
}


void Foam::solvers::multiPhysicsSolver::correctPhysics()
{
    scalar iterN(0);
    scalar residual(0);

    const Time& runTime(meshHandler_->returnMesh(singlePhysicsSolverNames_[0]).time());

    do
    {
        meshHandler_->mapTheseFields(runTime, singlePhysicsSolverNames_ );
        forAll(solvers_, solvI)
        {
            solvers_[solvI].deformMesh();
            (iterN>0 and correctOnlyEnergy_[solvI])?
                solvers_[solvI].correctTightlyCoupledPhysics():
                solvers_[solvI].correctPhysics();
            solvers_[solvI].correctBaffleLessFields();
        }

        forAll(solvers_, solvI)
        {
            residual = max(residual, solvers_[solvI].getResidual());
        }

        ++iterN;

        if(residual<minResidual_)
            Info << nl<<"Multiphysics loop converged after " << iterN <<" iterations"<<endl<<nl;
    }
    while(residual>minResidual_ && iterN < maxIterations_);
}

void Foam::solvers::multiPhysicsSolver::correctTightlyCoupledPhysics()
{
    correctPhysics();
}

scalar Foam::solvers::multiPhysicsSolver::getResidual()
{
    scalar residual = 0;
    forAll(solvers_, solvI)
    {
        residual = max(residual, solvers_[solvI].getResidual());
    }

    return residual;
}

scalar Foam::solvers::multiPhysicsSolver::maxDeltaT()
{
    scalar maxDeltaT= VGREAT;
    forAll(solvers_, solvI)
    {
        maxDeltaT = min(maxDeltaT, solvers_[solvI].maxDeltaT());
    }

    return maxDeltaT;
}


// ************************************************************************* //
