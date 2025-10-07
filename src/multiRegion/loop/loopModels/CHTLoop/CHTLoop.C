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

#include "CHTLoop.H"
#include "zeroGradientFvPatchFields.H"
#include "mergeOrSplitBaffles.H"
#include "fvmSup.H"
#include "fixedValuePointPatchFields.H"
#include "primitivePatchInterpolation.H"
#include "mappedPatchBase.H"
#include "mixedFvPatchField.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(CHTLoop, 0);
    addToRunTimeSelectionTable(solver, CHTLoop, dynamicFvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::CHTLoop::CHTLoop
(
    dynamicFvMesh& mesh 
)
:
    loop(mesh)
{
}


// * * * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //
void Foam::solvers::CHTLoop::createSolvers(word name)
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

    CHTProperties_ = multiPhysicsDict_.subDict("CHTLoopCoeffs");
    fluidRegionName_ = CHTProperties_.get<word>("fluidRegionName");
    solidRegionName_ = CHTProperties_.get<word>("solidRegionName");
    couplingStartTime_ = CHTProperties_.getOrDefault<scalar>("couplingStartTime",0);
    minResidual_ = multiPhysicsDict_.get<scalar>("minResidual");
    maxIterations_ = multiPhysicsDict_.get<label>("maxIterations");
    useHTC_ = CHTProperties_.getOrDefault<bool>("useHTC", false);
    oneWayCoupling_ = CHTProperties_.getOrDefault<bool>("oneWayCoupling", false);
    fluidKappa_ = CHTProperties_.getOrDefault<word>("fluidKappa", "kappaEff");
    solidKappa_ = CHTProperties_.getOrDefault<word>("solidKappa", "k");


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
    

    // Get list of patches 
    wordList fluidPatches = CHTProperties_.get<wordList>("fluidPatches");
    wordList solidPatches = CHTProperties_.get<wordList>("solidPatches");
    fluidPatchIDs_.setSize(fluidPatches.size());
    solidPatchIDs_.setSize(solidPatches.size());

    if(fluidPatches.size() != solidPatches.size())
    {
        FatalErrorInFunction
        << "Number of fluid and solid patches must be equal"
        <<  endl
        << "Fluid patches: " << fluidPatches.size() << endl
        << "Solid patches: " << solidPatches.size() << endl
        << exit(FatalError);
    }

    forAll(fluidPatches, wordI)
    {
        label fluidPatchID = meshHandler_->returnMesh(fluidRegionName_).boundaryMesh().findPatchID(fluidPatches[wordI]);
        if(fluidPatchID == -1)
        {
            FatalErrorInFunction
            << "Patch " << fluidPatches[wordI] << " not found in fluid region "
            << fluidRegionName_ << endl
            << exit(FatalError);
        }
        fluidPatchIDs_[wordI]= fluidPatchID;

        label solidPatchID = meshHandler_->returnMesh(solidRegionName_).boundaryMesh().findPatchID(solidPatches[wordI]);
        if(solidPatchID == -1)
        {
            FatalErrorInFunction
            << "Patch " << fluidPatches[wordI] << " not found in solid region "
            << solidRegionName_ << endl
            << exit(FatalError);
        }
        solidPatchIDs_[wordI] = solidPatchID;
    }
    
    
}

void Foam::solvers::CHTLoop::correctPhysics()
{
    if(oneWayCoupling_)
    {
        word targetRegion  = CHTProperties_.get<word>("targetRegion");
        if(targetRegion == fluidRegionName_)
        {
            solvers_[solidSolverID_].correctPhysics();
            if(meshHandler_->returnMesh(solidRegionName_).time().value()>couplingStartTime_)
                fromSolidToFluid();
            solvers_[fluidSolverID_].correctPhysics();
        }
        else if(targetRegion == solidRegionName_)
        {
            solvers_[fluidSolverID_].correctPhysics();
            if(meshHandler_->returnMesh(solidRegionName_).time().value()>couplingStartTime_)
                fromFluidToSolid();
            solvers_[solidSolverID_].correctPhysics();            
        }
        else
        {
            FatalErrorInFunction
            << "Target region must be either "
            <<  fluidRegionName_
            << " or " 
            << solidRegionName_ << endl
            << exit(FatalError);
        }
    }
    
    else
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
}


void Foam::solvers::CHTLoop::fromFluidToSolid()
{

    dynamicFvMesh& solidMesh = const_cast<dynamicFvMesh&>(meshHandler_->returnMesh(solidRegionName_));

    // Get patch and field

    const volScalarField& Ts = solidMesh.lookupObject<volScalarField>("T");

    forAll(solidPatchIDs_, solidPatchI)
    {
        const fvPatchScalarField& T_patch = Ts.boundaryField()[solidPatchIDs_[solidPatchI]];

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

        bool useHTC = CHTProperties_.get<bool>("useHTC");

        scalarList htcFluid = 
        (
            useHTC ?
            nbrPatch.lookupPatchField<volScalarField, scalar>("htc") :
            nbrPatch.lookupPatchField<volScalarField , scalar>(fluidKappa_) * nbrPatch.deltaCoeffs()
        );

        scalarField solidWeight = 
        (
            solidMesh.boundary()[solidPatchIDs_[solidPatchI]].lookupPatchField<volScalarField, scalar>(solidKappa_)
            * solidMesh.boundary()[solidPatchIDs_[solidPatchI]].deltaCoeffs()
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

}

void Foam::solvers::CHTLoop::fromSolidToFluid()
{

    dynamicFvMesh& fluidMesh = const_cast<dynamicFvMesh&>(meshHandler_->returnMesh(fluidRegionName_));


    const volScalarField& Tf = fluidMesh.lookupObject<volScalarField>("T");

    forAll(fluidPatchIDs_, fluidPatchI)
    {
    
        const fvPatchScalarField& T_patch = Tf.boundaryField()[fluidPatchIDs_[fluidPatchI]];

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

        scalarList htcSolid = 
            nbrPatch.lookupPatchField<volScalarField , scalar>(solidKappa_) * nbrPatch.deltaCoeffs();
        
        bool useHTC(CHTProperties_.get<bool>("useHTC"));
        scalarField fluidWeight = 
        (
            useHTC ?
            (fluidMesh.boundary()[fluidPatchIDs_[fluidPatchI]].lookupPatchField<volScalarField, scalar>("htc")) :
            (fluidMesh.boundary()[fluidPatchIDs_[fluidPatchI]].lookupPatchField<volScalarField, scalar>(fluidKappa_)
            * fluidMesh.boundary()[fluidPatchIDs_[fluidPatchI]].deltaCoeffs())
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


}

 Foam::scalar Foam::solvers::CHTLoop::calcFSIResidual()
{
    return scalar(0);
}

// ************************************************************************* //
