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

#include "FSILoop.H"
#include "zeroGradientFvPatchFields.H"
#include "mergeOrSplitBaffles.H"
#include "fvmSup.H"
#include "fixedValuePointPatchFields.H"
#include "primitivePatchInterpolation.H"
#include "mappedPatchBase.H"
#include "tractionDisplacementFvPatchVectorField.H"
#include "mixedFvPatchField.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(FSILoop, 0);
    addToRunTimeSelectionTable(solver, FSILoop, dynamicFvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::FSILoop::FSILoop
(
    dynamicFvMesh& mesh 
)
:
    loop(mesh)
{
}


// * * * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //
void Foam::solvers::FSILoop::createSolvers(word name)
{
    multiPhysicsDict_ = this->subDict("regionSolvers").subDict(name);

    solverNames_ = multiPhysicsDict_.subDict("subSolvers").toc();
    singlePhysicsSolverNames_.setSize(0);

    // First create the new solvers

    solvers_.setSize(solverNames_.size());

    wordList energyOnlyList = multiPhysicsDict_.getOrDefault<wordList>("includeFluidMechanicsInLoop", wordList());

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

        static const HashSet<word> loopTypes = { "picardLoop", "picardLoopNoFluid", "FSILoop", "CHTLoop", "multiScaleLoop"};

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

    FSIProperties_ = multiPhysicsDict_.subDict("FSILoopCoeffs");
    fluidRegionName_ = FSIProperties_.get<word>("fluidRegionName");
    solidRegionName_ = FSIProperties_.get<word>("solidRegionName");
    underRelaxation_ = FSIProperties_.get<scalar>("underRelaxationFactor");
    couplingStartTime_ = FSIProperties_.getOrDefault<scalar>("couplingStartTime",0);
    thermalCoupling_ = FSIProperties_.getOrDefault<bool>("thermalCoupling",false);
    minResidual_ = multiPhysicsDict_.get<scalar>("minResidual");
    maxIterations_ = multiPhysicsDict_.get<label>("maxIterations");


        // Get indices of solvers

    if(solvers_.size()!=2)
    {
        FatalErrorInFunction
        << "Please select one fluid and one solid solver"
        <<  endl
        << "Valid types are: "  << endl
        <<"Fluid: 2(rhoPimpleFoam, onePhase) " << endl
        <<"Solid: 2(extendedThermoMechanics, legacyThermoMechanics)" << endl
        << exit(FatalError);
    }
    else
    {
        if(solvers_[0].regionName() == fluidRegionName_)
        {
            fluidSolverID_ = 0;
        }
        else
        {
            fluidSolverID_ = 1;
        }
        solidSolverID_ = 1-fluidSolverID_;

        word fluidSolverType = 
            multiPhysicsDict_.subDict("subSolvers").
            get<word>(fluidRegionName_);
        
        word solidSolverType = 
            multiPhysicsDict_.subDict("subSolvers").
            get<word>(solidRegionName_);

        if 
        (
            (fluidSolverType != "onePhase" and fluidSolverType != "rhoPimpleFoam") 
            or
            (solidSolverType != "legacyThermomechanics" and solidSolverType != "extendedThermoMechanics" and solidSolverType != "fuelBehaviour")

        )
        {
            FatalErrorInFunction
            << "Please select one fluid and one solid solver"
            <<  endl
            << "Valid types are: "  << endl
            <<"Fluid: 2(rhoPimpleFoam, onePhase) " << endl
            <<"Solid: 3(extendedThermoMechanics, legacyThermoMechanics, fuelBehaviour)" << endl
            << exit(FatalError);
        }

    }
    
    word fluidSidePatchName = FSIProperties_.get<word>("fluidSidePatchName");
    fluidPatchID_ = meshHandler_->returnMesh(fluidRegionName_).boundaryMesh().findPatchID(fluidSidePatchName);
    word solidSidePatchName = FSIProperties_.get<word>("solidSidePatchName");
    solidPatchID_ = meshHandler_->returnMesh(solidRegionName_).boundaryMesh().findPatchID(solidSidePatchName);


    solidInterpolator_.reset
    (
        new volPointInterpolation
        (
            meshHandler_->returnMesh(solidRegionName_)
        )
    );



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

void Foam::solvers::FSILoop::correctPhysics()
{
    // Start loop 

    scalar iterN(0);
    scalar residual(0);

    do 
    {
        // fromSolidToFluid();
        solvers_[fluidSolverID_].correctPhysics();
        if(meshHandler_->returnMesh(solidRegionName_).time().value()>couplingStartTime_)
            fromFluidToSolid();
        solvers_[solidSolverID_].correctPhysics();
        if(meshHandler_->returnMesh(solidRegionName_).time().value()>couplingStartTime_)
            fromSolidToFluid();

        
        if(meshHandler_->returnMesh(solidRegionName_).time().value()>2)
            residual = calcFSIResidual();

        ++iterN;

        if(residual<minResidual_)
            Info << nl<<"Multiphysics loop converged after " << iterN <<" iterations"<<endl<<nl;
    }
    while(residual>minResidual_ && iterN < maxIterations_); 
}


void Foam::solvers::FSILoop::fromFluidToSolid()
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


    // Optional step - if thermal coupling map T and q
    if(thermalCoupling_)
    {

        // Get patch and field

        const volScalarField& Ts = solidMesh.lookupObject<volScalarField>("T");

        const fvPatchScalarField& T_patch = Ts.boundaryField()[solidPatchID_];

        const mixedFvPatchField<scalar>& mixedPatch =
        refCast<const mixedFvPatchField<scalar>>(T_patch);

        scalarField& refVal = const_cast<scalarField&>(mixedPatch.refValue());
        scalarField& valueFrac = const_cast<scalarField&>(mixedPatch.valueFraction());
        scalarField& refGrad = const_cast<scalarField&>(mixedPatch.refGrad());


        int oldTag = UPstream::msgType();
        UPstream::msgType() = oldTag+1;


        const mappedPatchBase& Tmpp = refCast<const mappedPatchBase>
        (
            mixedPatch.patch().patch()
        );

        const polyMesh& nbrMesh = Tmpp.sampleMesh();
        const fvPatch& nbrPatch = refCast<const fvMesh>
        (
            nbrMesh
        ).boundary()[Tmpp.samplePolyPatch().index()];

        scalarList temperatureFluid =
        nbrPatch.lookupPatchField<volScalarField, scalar>("T");

        bool useHTC = FSIProperties_.get<bool>("useHTC");
        word fluidKappa = FSIProperties_.getOrDefault<word>("fluidKappa", "kappaEff");
        word solidKappa = FSIProperties_.getOrDefault<word>("solidKappa", "k");

        scalarList htcFluid = 
        (
            useHTC ?
            nbrPatch.lookupPatchField<volScalarField, scalar>("htc") :
            nbrPatch.lookupPatchField<volScalarField , scalar>(fluidKappa) * nbrPatch.deltaCoeffs()
        );

        scalarField solidWeight = 
        (
            solidMesh.boundary()[solidPatchID_].lookupPatchField<volScalarField, scalar>(solidKappa)
            * solidMesh.boundary()[solidPatchID_].deltaCoeffs()
        );

        // Set ref value

        Tmpp.distribute(temperatureFluid);
        refVal = temperatureFluid;

        //Set value fraction
        Tmpp.distribute(htcFluid);
        valueFrac = htcFluid / (solidWeight + htcFluid);

        //Set gradient
        refGrad = 0;

    }

    // Step 3 - Deform back the solid mesh


    displacementPoints -= pointD;

    solidMesh.movePoints(displacementPoints);


}

void Foam::solvers::FSILoop::fromSolidToFluid()
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

    // Optional step - if thermal coupling map T and q
    if(thermalCoupling_)
    {

        // Get patch and field

        const volScalarField& Tf = fluidMesh.lookupObject<volScalarField>("T");

        const fvPatchScalarField& T_patch = Tf.boundaryField()[fluidPatchID_];

        const mixedFvPatchField<scalar>& mixedPatch =
        refCast<const mixedFvPatchField<scalar>>(T_patch);

        scalarField& refVal = const_cast<scalarField&>(mixedPatch.refValue());
        scalarField& valueFrac = const_cast<scalarField&>(mixedPatch.valueFraction());
        scalarField& refGrad = const_cast<scalarField&>(mixedPatch.refGrad());


        int oldTag = UPstream::msgType();
        UPstream::msgType() = oldTag+1;


        const mappedPatchBase& Tmpp = refCast<const mappedPatchBase>
        (
            mixedPatch.patch().patch()
        );

        const polyMesh& nbrMesh = Tmpp.sampleMesh();
        const fvPatch& nbrPatch = refCast<const fvMesh>
        (
            nbrMesh
        ).boundary()[Tmpp.samplePolyPatch().index()];

        scalarList temperatureSolid =
        nbrPatch.lookupPatchField<volScalarField, scalar>("T");

        bool useHTC(FSIProperties_.get<bool>("useHTC"));
        word fluidKappa = FSIProperties_.getOrDefault<word>("fluidKappa", "kappaEff");
        word solidKappa = FSIProperties_.getOrDefault<word>("solidKappa", "k");

        scalarList htcSolid = 
            nbrPatch.lookupPatchField<volScalarField , scalar>(solidKappa) * nbrPatch.deltaCoeffs();

        scalarField fluidWeight = 
        (
            useHTC ?
            (fluidMesh.boundary()[fluidPatchID_].lookupPatchField<volScalarField, scalar>("htc")) :
            (fluidMesh.boundary()[fluidPatchID_].lookupPatchField<volScalarField, scalar>(fluidKappa)
            * fluidMesh.boundary()[fluidPatchID_].deltaCoeffs())
        );

        // Set ref value

        Tmpp.distribute(temperatureSolid);
        refVal = temperatureSolid;

        //Set value fraction
        Tmpp.distribute(htcSolid);
        valueFrac = htcSolid / (fluidWeight + htcSolid);

        //Set gradient
        refGrad = 0;

    }

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

 Foam::scalar Foam::solvers::FSILoop::calcFSIResidual()
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

// ************************************************************************* //
