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
            runTime.time().system(),
            runTime.db(),
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

/*void Foam::solvers::loop::PicardIterRoutine()
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

void Foam::solvers::loop::FSIInitialization()
{
    underRelaxation_ = multiPhysicsDict_.get<scalar>("underRelaxation");
    fluidRegionName_ = multiPhysicsDict_.get<word>("fluidRegionName");
    solidRegionName_ = multiPhysicsDict_.get<word>("solidRegionName");
    word fluidSidePatchName = multiPhysicsDict_.get<word>("fluidSidePatchName");
    fluidPatchID_ = meshHandler_->returnMesh(fluidRegionName_).boundaryMesh().findPatchID(fluidSidePatchName);
    word solidSidePatchName = multiPhysicsDict_.get<word>("solidSidePatchName");
    solidPatchID_ = meshHandler_->returnMesh(solidRegionName_).boundaryMesh().findPatchID(solidSidePatchName);


    // Get patch interpolator

    primitivePatchInterpolation patchInterpolator
    (
        meshHandler_->returnMesh(fluidRegionName_).boundaryMesh()[fluidPatchID_]
    );

    // Get solid displacement field

    vectorField displacementAtFaces
    (
        meshHandler_->returnMesh(solidRegionName_).lookupObject<volVectorField>("D").
        boundaryField()[solidPatchID_]
    );


    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
        meshHandler_->returnMesh(fluidRegionName_).boundaryMesh()[fluidPatchID_]
    );

    mpp.distribute(displacementAtFaces);


    oldDisplacementAtFaces_.reset
    (
        new vectorField
        (
            displacementAtFaces
        )
    );

    if(solverNames_.size()>2)
    {
        FatalErrorInFunction << "Please select max 2 solvers" << endl;

    }

    forAll(solverNames_,nameI)
    {
        word solverType
        (
            multiPhysicsDict_
                .subDict("subSolvers")
                .get<word>(solverNames_[nameI])
        );

        if(solverType == "pimpleFluid")
        {
            fluidSolverID_=nameI;
        }
        else if (solverType == "extendedThermoMechanics")
        {
            solidSolverID_=nameI;
        }
        else
        {
            FatalErrorInFunction << "Please select FSI solvers" << endl << exit(FatalError);
        }
    }

    solidInterpolator_.reset
    (
        new volPointInterpolation
        (
            meshHandler_->returnMesh(solidRegionName_)
        )
    );

    fluidDisplacementAtSolid_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "fluidDisplacementAtSolid",
                meshHandler_->returnMesh(solidRegionName_).time().timeName(),
                meshHandler_->returnMesh(solidRegionName_),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            meshHandler_->returnMesh(solidRegionName_),
            dimensionedVector("fluidDisplacementAtSolid", dimLength, vector::zero)
        )
    );

}


void Foam::solvers::loop::FSIRoutine()
{
    // To do's:
        // Get list of coupled patches - OK
        // Make sure that the patches are of the right type (i.e., movingWallVelocity, fixedGradient and mappedTraction)
        // All the fields mappings are handled at the patch level so not much to do there - OK
        // Take care of fluid deformation (this also implicates using dynamicFvMesh (OK) and storing displacement fields (OK)) OK


    // Start loop

    scalar iterN(0);
    scalar residual(0);

    do
    {
        // fromSolidToFluid();
        solvers_[fluidSolverID_].correctPhysics();
        if(meshHandler_->returnMesh(solidRegionName_).time().value()>2)
            fromFluidToSolid();
        solvers_[solidSolverID_].correctPhysics();
        if(meshHandler_->returnMesh(solidRegionName_).time().value()>2)
            fromSolidToFluid();


        if(meshHandler_->returnMesh(solidRegionName_).time().value()>2)
            residual = calcFSIResidual();

        ++iterN;

        if(residual<minResidual_)
            Info << nl<<"Multiphysics loop converged after " << iterN <<" iterations"<<endl<<nl;
    }
    while(residual>minResidual_ && iterN < maxIterations_);

}

// void Foam::solvers::loop::deformFluidMeshForFSI()
// {
//     dynamicFvMesh& fluidMesh = const_cast<dynamicFvMesh&>(meshHandler_->returnMesh(fluidRegionName_));
//     dynamicFvMesh& solidMesh = const_cast<dynamicFvMesh&>(meshHandler_->returnMesh(solidRegionName_));

//     pointVectorField& motionU = const_cast<pointVectorField&>
//     (
//         fluidMesh.objectRegistry::
//         lookupObject<pointVectorField>
//         (
//             "pointMotionU"
//         )
//     );

//     fixedValuePointPatchVectorField& motionUFluidPatch = refCast<fixedValuePointPatchVectorField>
//     (
//         motionU.boundaryFieldRef()[fluidPatchID_]
//     );

//     // Get patch interpolator

//     primitivePatchInterpolation patchInterpolator
//     (
//         fluidMesh.boundaryMesh()[fluidPatchID_]
//     );

//     // Get solid displacement field

//     vectorField displacementAtFaces
//     (
//         solidMesh.lookupObject<volVectorField>("D").
//         boundaryField()[solidPatchID_]
//     );

//     // Map solid displacement to fluid mesh

//     int oldTag = UPstream::msgType();
//     UPstream::msgType() = oldTag+1;

//     // Get the coupling information from the mappedPatchBase
//     const mappedPatchBase& mpp = refCast<const mappedPatchBase>
//     (
//         fluidMesh.boundaryMesh()[fluidPatchID_]
//     );

//     mpp.distribute(displacementAtFaces);

//     vectorField deltaDAtPoints =
//     patchInterpolator.faceToPointInterpolate
//     (
//         displacementAtFaces - oldDisplacementAtFaces_()
//     );


//     motionUFluidPatch == deltaDAtPoints /fluidMesh.time().deltaT().value();

//     fluidMesh.update();
// }

void Foam::solvers::loop::fromFluidToSolid()
{
    // Steps to do:
        //Deform the solid mesh
        //Pass field
        //Undeform solid mesh

    // Step 1 - Deform solid mesh

    dynamicFvMesh& solidMesh = const_cast<dynamicFvMesh&>(meshHandler_->returnMesh(solidRegionName_));

    pointField displacementPoints = solidMesh.points();

    tmp<pointVectorField> pointD
    (
        solidInterpolator_().interpolate
        (
            fluidDisplacementAtSolid_()
        )
    );


    displacementPoints += pointD;

    solidMesh.movePoints(displacementPoints);

    // Step 2 - Pass info from fluid to solid

    const volVectorField& D = solidMesh.lookupObject<volVectorField>("D");

    const fvPatchVectorField& D_patch = D.boundaryField()[solidPatchID_];

    const tractionDisplacementFvPatchVectorField& tractionPatch =
    refCast<const tractionDisplacementFvPatchVectorField>(D_patch);

    vectorField& traction = const_cast<vectorField&>(tractionPatch.traction());

    // Get the coupling information from the mappedPatchBase


    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;


    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
        tractionPatch.patch().patch()
    );

    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch = refCast<const fvMesh>
    (
        nbrMesh
    ).boundary()[mpp.samplePolyPatch().index()];



    scalarList pressureFluid =
        nbrPatch.lookupPatchField<volScalarField, scalar>("p");

    //mpp.distribute(pressureFluid);

    vectorField tractionFluid =
    (
        (nbrPatch.lookupPatchField<volSymmTensorField, symmTensor>("stressTensor")
        & nbrPatch.nf())
    );

    mpp.distribute(tractionFluid);


    //pressure = pressureFluid;
    traction = -tractionFluid ;

    // Step 3 - Deform back the solid mesh


    displacementPoints -= pointD;

    solidMesh.movePoints(displacementPoints);


}

void Foam::solvers::loop::fromSolidToFluid()
{

    // Steps to do:
    // Deform the solid mesh with old displacement
    // Pass field
    // Undeform solid mesh with old displacement
    // Deform fluid mesh
    // Update displacement



    // Step 1 - Deform solid mesh with old displacement


    dynamicFvMesh& solidMesh = const_cast<dynamicFvMesh&>(meshHandler_->returnMesh(solidRegionName_));

    pointField displacementPoints = solidMesh.points();

    const volVectorField& solidDisplacement = solidMesh.lookupObject<volVectorField>("D");

    tmp<pointVectorField> pointD
    (
        solidInterpolator_().interpolate
        (
            fluidDisplacementAtSolid_()
        )
    );

    displacementPoints += pointD;

    solidMesh.movePoints(displacementPoints);

    // Step 2 - Map displacement from solid to fluid

    dynamicFvMesh& fluidMesh = const_cast<dynamicFvMesh&>(meshHandler_->returnMesh(fluidRegionName_));

    vectorField displacementAtFaces
    (
        solidMesh.lookupObject<volVectorField>("D").
        boundaryField()[solidPatchID_]
    );

    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
        fluidMesh.boundaryMesh()[fluidPatchID_]
    );

    mpp.distribute(displacementAtFaces);

    // Step 3 - Un-deform solid mesh


    displacementPoints -= pointD;

    solidMesh.movePoints(displacementPoints);

    // Step 4 - Update the dynamicFvMesh solver

    pointVectorField& motionU = const_cast<pointVectorField&>
    (
        fluidMesh.objectRegistry::
        lookupObject<pointVectorField>
        (
            "pointMotionU"
        )
    );

    fixedValuePointPatchVectorField& motionUFluidPatch = refCast<fixedValuePointPatchVectorField>
    (
        motionU.boundaryFieldRef()[fluidPatchID_]
    );

    // Get patch interpolator

    primitivePatchInterpolation patchInterpolator
    (
        fluidMesh.boundaryMesh()[fluidPatchID_]
    );

    displacementAtFaces = displacementAtFaces *underRelaxation_ + oldDisplacementAtFaces_()*(1-underRelaxation_);

    vectorField deltaDAtPoints =
    patchInterpolator.faceToPointInterpolate
    (
        displacementAtFaces - oldDisplacementAtFaces_()
    );


    motionUFluidPatch == deltaDAtPoints /fluidMesh.time().deltaT().value();

    fluidMesh.update();

    // Step 5 - Update displacements

    fluidDisplacementAtSolid_() = solidDisplacement*underRelaxation_ + (1-underRelaxation_)*fluidDisplacementAtSolid_();

    oldDisplacementAtFaces_() = displacementAtFaces;
}

 Foam::scalar Foam::solvers::loop::calcFSIResidual()
{
    // Steps to do:
        //Step 1: Deform solid mesh
        //Step 2: Map solid displacement onto fluid mesh
        //Step 3: Compare the solid to the fluid displacement
        //Step 4: Compute the residual
        //Step 5: Un-deform solid

    //Step 1: Deform solid

    dynamicFvMesh& solidMesh = const_cast<dynamicFvMesh&>(meshHandler_->returnMesh(solidRegionName_));

    pointField displacementPoints = solidMesh.points();

    tmp<pointVectorField> pointD
    (
        solidInterpolator_().interpolate
        (
            fluidDisplacementAtSolid_()
        )
    );

    displacementPoints += pointD;

    solidMesh.movePoints(displacementPoints);

    // Step 2 - Map displacement from solid to fluid

    dynamicFvMesh& fluidMesh = const_cast<dynamicFvMesh&>(meshHandler_->returnMesh(fluidRegionName_));

    const volVectorField& D = solidMesh.lookupObject<volVectorField>("D");

    vectorField displacementAtFaces
    (
        D.
        boundaryField()[solidPatchID_]
    );

    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
        fluidMesh.boundaryMesh()[fluidPatchID_]
    );

    mpp.distribute(displacementAtFaces);

    // Step 3 - Un-deform solid mesh


    displacementPoints -= pointD;

    solidMesh.movePoints(displacementPoints);

    // Step 4 - Compare fluid displacement to solid

    vectorField residual = displacementAtFaces-oldDisplacementAtFaces_();

    Info << "New residual: " << Foam::sqrt(gSum(magSqr(residual))) << endl;

    //Return normalised residual
    return
    (
        Foam::sqrt(gSum(magSqr(residual)))
    );

}

*/


// ************************************************************************* //
