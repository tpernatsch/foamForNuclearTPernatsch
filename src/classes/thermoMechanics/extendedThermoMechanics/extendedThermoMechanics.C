/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "extendedThermoMechanics.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

    defineTypeNameAndDebug(extendedThermoMechanics, 0);
    addToRunTimeSelectionTable
    (
        solver, 
        extendedThermoMechanics, 
        dynamicFvMesh
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::solvers::extendedThermoMechanics::extendedThermoMechanics
(
    dynamicFvMesh& mesh
)
:
    thermoMechanics(mesh),
    runTime_(static_cast<const offbeatTime&>(mesh.time())),
    mesh_(mesh),
    solverDict_
    (
        IOobject
        (
            "thermoMechanicalProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    adjustableTime_(mesh_, const_cast<Time&>(mesh_.time())),
    globalOpt_(mesh_, solverDict_),
    mat_(materials::New(mesh_, solverDict_)),
    mapper_(sliceMapper::New(mesh_, mat_, solverDict_)),
    rheo_(rheology::New(mesh_, mat_, solverDict_)),
    heatSrc_(heatSource::New(mesh_, mat_, solverDict_)),
    thermal_(thermalSubSolver::New(mesh_, mat_(), solverDict_)),
    mechanics_(mechanicsSubSolver::New(mesh_, mat_, rheo_(), solverDict_))
{

    // Check user parameters
    listUserParameters();
    listRegisteredUserParameters();
    checkUserParameters();

    // Set last time marker to adjustable time step
    adjustableTime_.setLastTimeMarker
    (
        min
        (
            heatSrc_->lastTimeMarker(),
            runTime_.endTime().value()
        )
    );

}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::solvers::extendedThermoMechanics::correctPhysics()
{

    mesh_.update();

    Info << "Solving region: " << mesh_.name()<<nl<<endl;
    // True when time step is converged
    bool converged_ = false;

    // Read from fvSolution dict the maximum number of outer iterations
    const dictionary& stressControl = 
    mesh_.solutionDict().subDict("stressAnalysis");
    int maxOuterIter(readInt(stressControl.lookup("maxOuterIter")));

    // Track number of outer iterations (per time-step)
    int nOuterIter = 0;

    // Update mesh if necessary
    if
    (
        mesh_.foundObject<fvMesh>("referenceMesh")
    )
    {
        mechanics_->updateMesh();
    }

    // Update incremental fields (DD and gradDD = 0) if necessary
    mechanics_->updateIncrementalFields();

    do 
    {
        Info << "OuterIteration n. " << nOuterIter << nl <<  endl;

        storeGlobalFieldsPrevIter();

        // First update corrosion, only then update AMI. Necessary because
        // if corrosion moves the mesh, the AMI is also updated given that 
        // moving the mesh clears the geometry and the AMIPtr. However this
        // first update, does not use the updated mesh location
        // TODO: find a way to update mesh without AMI?

        forAll(mesh_.boundary(), patchi)
        {
            if(isType<regionCoupledOFFBEATFvPatch>(mesh_.boundary()[patchi]))
            {
                const regionCoupledOFFBEATFvPatch& patch
                = refCast<const regionCoupledOFFBEATFvPatch>(mesh_.boundary()[patchi]);

                patch.updateAMI();
            }
        }

        bool isCouplingDict(solverDict_.found("couplingOptions"));

        bool getTFromTH
        (
            isCouplingDict ?
            solverDict_.subDict("couplingOptions").getOrDefault<bool>("correctTFromTH", false)
            : false
        );
        bool computeNeutroMeshDisplacement
        (
            isCouplingDict ?
            solverDict_.subDict("couplingOptions").getOrDefault<bool>("correctDispForNeutro", false)
            : false
        );

        heatSrc_->correct();
        thermal_->correct();
        if (getTFromTH)
        {
            #include "correctT.H"
        }
        mat_->correctBehavioralModels();
        rheo_->getAdditionalStrain();
        mechanics_->correct();
        if (computeNeutroMeshDisplacement)
        {
            #include "correctMeshDisp.H"
        }
        correctBCsGlobalFields();
        relaxGlobalFields();

        nOuterIter++;
        
        // Check convergence
        converged_ = 
        (
            mechanics_->converged() 
            and 
            thermal_->converged()
            and 
            rheo_->converged()
        );
    }
    while
    (
        (nOuterIter < 2)  
        || 
        (not(converged_) && nOuterIter < maxOuterIter) 
    );

    mechanics_->updateTotalFields();

    // Check if failure occurred
    mat_->checkFailure();

}

void Foam::solvers::extendedThermoMechanics::correctTightlyCoupledPhysics()
{
    correctPhysics();
}

Foam::scalar Foam::solvers::extendedThermoMechanics::maxDeltaT()
{   
    scalar deltaT(GREAT);
    deltaT = min(deltaT, mechanics_->nextDeltaT());
    deltaT = min(deltaT, heatSrc_->nextDeltaT());
    deltaT = min(deltaT, mat_->nextDeltaT());
    
    return deltaT;
}

Foam::scalar Foam::solvers::extendedThermoMechanics::getMaxToutClad
(
    const word& cladOuterPatchName
)
{
    const volScalarField Tref(mesh_.lookupObject<volScalarField>("T"));

    // Init Tmax
    scalar Tmax(0.0);

    forAll(mesh_.boundaryMesh(), patchI)
    {
        // Look for cladding outer patch 
        if( mesh_.boundaryMesh()[patchI].name() == cladOuterPatchName ) 
        {
            Tmax = max(Tref.boundaryField()[patchI]);
            break;
        }
    }

    return Tmax;
}


// ************************************************************************* //
