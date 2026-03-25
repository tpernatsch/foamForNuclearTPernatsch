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
#include "globalFieldLists.H"


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
    solver(mesh),
    runTime_(refCast<const offbeatTime>(mesh.time())),
    solverDict_
    (
        IOobject
        (
            "solverDict",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    globalOpt_(mesh, solverDict_),
    mat_(materials::New(mesh, solverDict_)),
    mapper_(sliceMapper::New(mesh, mat_(), solverDict_)),
    rheo_(rheology::New(mesh, mat_(), solverDict_)),
    heatSrc_(heatSource::New(mesh, mat_(), solverDict_)),
    burnup_(burnup::New(mesh, mat_(), solverDict_)),
    neutronics_(neutronicsSubSolver::New(mesh, burnup_, solverDict_)),
    fastFlux_(fastFlux::New(mesh, mat_(), solverDict_)),
    thermal_(thermalSubSolver::New(mesh, mat_(), solverDict_)),
    mechanics_(mechanicsSubSolver::New(mesh, mat_(), rheo_(), solverDict_)),
    fgr_(fissionGasRelease::New(mesh, mat_(), solverDict_)),    
    gapGas_(gapGasModel::New(mesh, mat_(), thermal_->T(), mechanics_->DorDD(), fgr_,  solverDict_)),
    corrosion_(corrosion::New(mesh, mat_(), solverDict_)),
    elementTransport_(elementTransport::New(mesh, mat_(), solverDict_)),
    #ifdef OPENFOAMFOUNDATION
    flow_(flowSubSolver::New(runTime_, solverDict_)),
    #endif
    nOuterIter_(0.0),
    totalIterations_(0.0)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
#ifdef isCommDataLayerIncluded
void Foam::solvers::extendedThermoMechanics::correctPhysics(fmuControl* fmu)
#else
void Foam::solvers::extendedThermoMechanics::correctPhysics()
#endif
{    
    #ifdef isCommDataLayerIncluded
    if (runTime_.timeIndex() == 0 && isSolveFMI)
    {
        fmu->receive();
        fmu->send();
    }
    #endif
    
    // True when time step is converged
    bool converged = false;

    // Read from fvSolution dict the maximum number of outer iterations
    const dictionary& stressControl = 
    mesh_.solutionDict().subDict("stressAnalysis");

    int maxOuterIter(readInt(stressControl.lookup("maxOuterIter")));
        
    int minOuterIter;
    if (maxOuterIter == 1)
    {
        // The user forces a single outer iteration per time-step
        minOuterIter = 1;
    }
    else
    {
        // Forcing at least 2 outer iteration to take into account non-linear feedback
        // before residuals are converged
        minOuterIter = 2;
    }

    #ifdef isCommDataLayerIncluded
    do
    {
        if (isSolveFMI) fmu->receive();
    #endif

    // Reset number of outer iterations (per time-step)
    nOuterIter_ = 0;

    // Update mesh if necessary
    if
    (
        mesh().foundObject<fvMesh>("referenceMesh")
    )
    {
        mechanics_->updateMesh();
    }

    // Update mesh due to transport solvers requirements if needed
    elementTransport_->updateMesh();
    
    // Update incremental fields (DD and gradDD = 0) if necessary
    mechanics_->updateIncrementalFields();
    
    #ifdef OPENFOAMFOUNDATION
    //- Initialise the flow solver before the solver loop
    flow_->init();
    #endif

    do 
    {
        Info << "OuterIteration n. " << nOuterIter_ << nl <<  endl;

        storeGlobalFieldsPrevIter();

        // First update corrosion, only then update AMI. Necessary because
        // if corrosion moves the mesh, the AMI is also updated given that 
        // moving the mesh clears the geometry and the AMIPtr. However this
        // first update, does not use the updated mesh location
        // TODO: find a way to update mesh without AMI?
        corrosion_->correct();

        // if(nOuterIter_ < 1)
        // {
            forAll(mesh().boundary(), patchi)
            {             
                if(isType<regionCoupledOFFBEATFvPatch>(mesh().boundary()[patchi]))
                {
                const regionCoupledOFFBEATFvPatch& patch
                    = refCast<const regionCoupledOFFBEATFvPatch>(mesh().boundary()[patchi]);

                    patch.updateAMI();
                }
            }
        // }

        fgr_->correct();
        gapGas_->correct();
        heatSrc_->correct();
        neutronics_->correct();
        burnup_->correct();
        fastFlux_->correct();
        thermal_->correct();
        mat_->correctBehavioralModels();
        rheo_->getAdditionalStrain();
        mechanics_->correct();
        elementTransport_->correct();

        #ifdef OPENFOAMFOUNDATION
        flow_->correct();
        #endif

        correctBCsGlobalFields();
        relaxGlobalFields();

        nOuterIter_++;
        totalIterations_++;
        
        // Check convergence
        converged = 
        (
            mechanics_->converged() 
            and 
            thermal_->converged()
            and 
            rheo_->converged()
            and
            neutronics_->converged()
            and
            elementTransport_->converged()
        #ifdef OPENFOAMFOUNDATION
            and
            flow_->converged()
        #endif
        );
    }
    while
    (
        (nOuterIter_ < minOuterIter)  
        || 
        (not(converged) && nOuterIter_ < maxOuterIter) 
    );

    #ifdef isCommDataLayerIncluded
    if (isSolveFMI)
    {
        commDataLayer& data = commDataLayer::New(runTime_);
        scalar& isConvergedFMI = data.getObj<scalar>
        (
            "isConverged",
            commDataLayer::causality::out
        );
        isConvergedFMI = converged ? 1.0 : 0.0;
    }

    if (isSolveFMI) fmu->send();
    } // End fmu implicit loop
    while (isSolveFMI && fmu->loop());
    #endif

    mechanics_->updateTotalFields();
    gapGas_->updateVariables();    
    fgr_->updateVariables();
    #ifdef OPENFOAMFOUNDATION
    flow_->finalise();
    #endif
    
    // Check if failure occurred
    mat_->checkFailure();
}


Foam::scalar Foam::solvers::extendedThermoMechanics::maxDeltaT()
{   
    scalar deltaT(GREAT);

    deltaT = min(deltaT, mechanics_->nextDeltaT());
    deltaT = min(deltaT, thermal_->nextDeltaT());
    deltaT = min(deltaT, heatSrc_->nextDeltaT());
    deltaT = min(deltaT, burnup_->nextDeltaT());
    deltaT = min(deltaT, fgr_->nextDeltaT());
    deltaT = min(deltaT, mat_->nextDeltaT());
    deltaT = min(deltaT, elementTransport_->nextDeltaT());
    #ifdef OPENFOAMFOUNDATION
    deltaT = min(deltaT, flow_->nextDeltaT());
    
    // Check for consistency of deltaT for parallel runs
    // If check fails then this suggests a coding error in the nextDeltaT
    // subroutine of one of the submodules
    if (Time::debug || offbeatTime::debug)
    {
        if (Pstream::parRun())
        {
            scalarField allDeltaT = helperFuns::gather(deltaT);

            if (Pstream::master())
            {
                if (mag(average(allDeltaT) - max(allDeltaT))/allDeltaT[0] > ROOTSMALL)
                {
                    FatalErrorInFunction()
                        << "Inconsistent deltaT... this suggests a coding error in the "
                        << "nextDeltaT subroutine of one of the submodules. Please submit "
                        << "a bug report." << nl
                        << tab << "Mechanics deltaT: " << helperFuns::gather(mechanics_->nextDeltaT()) << nl
                        << tab << "Thermal deltaT: " << helperFuns::gather(thermal_->nextDeltaT()) << nl
                        << tab << "Heat source deltaT: " << helperFuns::gather(heatSrc_->nextDeltaT()) << nl
                        << tab << "Burnup deltaT: " << helperFuns::gather(burnup_->nextDeltaT()) << nl
                        << tab << "FGR deltaT: " << helperFuns::gather(fgr_->nextDeltaT()) << nl
                        << tab << "Material deltaT: " << helperFuns::gather(mat_->nextDeltaT()) << nl
                        << tab << "Element transport deltaT: " << helperFuns::gather(elementTransport_->nextDeltaT()) << nl
                        << tab << "Flow deltaT: " << helperFuns::gather(flow_->nextDeltaT()) << nl
                        << tab << "Final deltaT: " << helperFuns::gather(deltaT) << nl
                        << tab << "Current time: " << helperFuns::gather(runTime_.value()) << endl
                        << abort(FatalError);
                }
            }
        }
    }
    #endif 
    
    return deltaT;
}


Foam::scalar Foam::solvers::extendedThermoMechanics::endTime()
{
    return heatSrc_->lastTimeMarker();
}


// ************************************************************************* //
