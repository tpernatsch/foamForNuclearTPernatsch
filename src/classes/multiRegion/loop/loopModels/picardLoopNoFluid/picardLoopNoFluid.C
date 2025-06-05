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

#include "picardLoopNoFluid.H"
#include "zeroGradientFvPatchFields.H"
#include "mergeOrSplitBaffles.H"
#include "fvmSup.H"
#include "fixedValuePointPatchFields.H"
#include "primitivePatchInterpolation.H"
#include "mappedPatchBase.H"
#include "tractionDisplacementFvPatchVectorField.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(picardLoopNoFluid, 0);
    addToRunTimeSelectionTable(solver, picardLoopNoFluid, dynamicFvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::picardLoopNoFluid::picardLoopNoFluid
(
    dynamicFvMesh& mesh
)
:
    loop(mesh)
{
}


// * * * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

void Foam::solvers::picardLoopNoFluid::createSolvers(word name)
{
    multiPhysicsDict_ = this->subDict("regionSolvers").subDict(name);

    solverNames_ = multiPhysicsDict_.subDict("subSolvers").toc();
    singlePhysicsSolverNames_.setSize(0);

    // First create the new solvers

    solvers_.setSize(solverNames_.size());

    forAll(solverNames_, nameI)
    {
        Info<< "Creating sub-scale solver " << solverNames_[nameI] << nl
            << endl;

        word solverType
        (
            multiPhysicsDict_
                .subDict("subSolvers")
                .get<word>(solverNames_[nameI])
        );

        // This replicates the structure of regionSolvers. This allows to create
        // subSolvers which are multiPhysicsSolvers themselves and so on,
        // creating a tree-like structure

        static const HashSet<word> loopTypes = { "picardLoopNoFluid", "FSILoop", "CHTLoop"};

        const word& meshName = solverNames_[nameI];
        word meshToUse = loopTypes.found(solverType) ? "dummy" : meshName;

        solvers_.set
        (
            nameI,
            solver::New(solverType, meshHandler_->returnMesh(meshToUse))
        );

        if (Foam::isA<Foam::solvers::loop>(solvers_[nameI]))
        {
            Foam::solvers::loop& loopRef =
                Foam::refCast<Foam::solvers::loop>(solvers_[nameI]);

            loopRef.getMapper(*meshHandler_);
            loopRef.createSolvers(solverNames_[nameI]);
        }
        else
        {
            singlePhysicsSolverNames_.append(solverNames_[nameI]);
        }
    }

    //loopProperties_ = multiPhysicsDict_.subDict("loopProperties");
    couplingStartTime_ = multiPhysicsDict_.getOrDefault<scalar>("couplingStartTime",0);
    minResidual_ = multiPhysicsDict_.get<scalar>("minResidual");
    maxIterations_ = multiPhysicsDict_.get<label>("maxIterations");
        
}

void Foam::solvers::picardLoopNoFluid::correctPhysics()
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
            solvers_[solvI].correctTightlyCoupledPhysics();
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



// ************************************************************************* //
