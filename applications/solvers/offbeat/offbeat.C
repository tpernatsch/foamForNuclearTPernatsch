/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Application
    offbeat

Description
    OFFBEAT (OpenFOAM Fuel BEhavior Analysis Tool) is a multi-dimensional fuel 
    performance code based on the open-source C++ library OpenFOAM. It can be
    used both for the analysis of standard 1.5D cases and for more complex 2D 
    and 3D local effects.

    TODO list:
    [-] Include fgr in gapGas->correctMass
    [-] Verify that 2 min outerIter are necessary
    [-] Description in file header
    [-] Eliminate writeStressField?

\*---------------------------------------------------------------------------*/

// Include standard classes
#include "fvCFD.H"
#include "SortableList.H"

// Include main OFFBEAT classes
#include "globalFieldLists.H"
#include "globalOptions.H"
#include "offbeatTime.H"
#include "thermalSubSolver.H"
#include "mechanicsSubSolver.H"
#include "gapGasModel.H"
#include "fissionGasRelease.H"
#include "burnup.H"
#include "neutronicsSubSolver.H"
#include "heatSource.H"
#include "fastFlux.H"
#include "timeHandling.H"
#include "sliceMapper.H"
#include "materials.H"
#include "rheology.H"
#include "corrosion.H"
#include "userParameters.H"
#include "elementTransport.H"
#include "dynamicFvMesh.H"

#if defined __has_include
#  if __has_include(<commDataLayer.H>)
#    include <commDataLayer.H>
#    define isCommDataLayerIncluded
#  endif
#endif

#ifdef isCommDataLayerIncluded
#include "commDataLayer.H"
#include "fmuControl.H"
#endif

#ifdef OPENFOAMFOUNDATION
#include "flowSubSolver.H"
#include "helperFuns.H"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Check if user requested post-processing 
int checkPostProcess(int argc, char *argv[])
{
    argList::addBoolOption
    (
        argList::postProcessOptionName,
        "Execute functionObjects only"
    );
    
    if (argList::postProcess(argc, argv))
    {
        Foam::timeSelector::addOptions();
        #include "addFunctionObjectOptions.H"

        // Set functionObject post-processing mode
        functionObject::postProcess = true;

        Foam::argList args(argc, argv);
        
        if (args.optionFound("list"))
        {
            functionObjectList::list();
            return 0;
        }
    }
    
    return 0;
}

//- If post-processing was requested execute all functionObjects
int postProcess(const argList& args, const Time& runTime)
{
    // Externally stored dictionary for functionObjectList
    // if not constructed from runTime
    dictionary controlDict = runTime.controlDict();

#ifdef OPENFOAMFOUNDATION
    HashSet<word> selectedFields;
#elif OPENFOAMESI
    HashSet<wordRe> selectedFields;
#endif

    // Construct functionObjectList
    autoPtr<functionObjectList> functionsPtr
    (
        functionObjectList::New
        (
            args,
            runTime,
            controlDict,
            selectedFields
        )
    );

    functionsPtr->execute();

    Info<< endl;

    Info<< "End\n" << endl;

    return 0;
}

int main(int argc, char *argv[])
{
    checkPostProcess(argc, argv);
    
#   include "setRootCase.H"

    Info<< "Create time\n" << Foam::endl;
  
    offbeatTime runTime
    (
        Foam::Time::controlDictName, 
        args
    );
 
#   include "createDynamicFvMesh.H"

    // Read main OFFBEAT dictionary
    IOdictionary solverDict
    (
        IOobject
        (
            "solverDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Create objects and autoPtr to main OFFBEAT classes
    timeHandling adjustableTime(mesh, runTime);
    globalOptions globalOpt_(mesh, solverDict);
    autoPtr<materials> mat_(materials::New(mesh, solverDict));
    materials& mat = mat_();
    
    autoPtr<sliceMapper> mapper_(sliceMapper::New(mesh, mat, solverDict));

    autoPtr<rheology> rheo_(rheology::New(mesh, mat, solverDict));    
    autoPtr<heatSource> heatSrc_(heatSource::New(mesh, mat, solverDict));
    autoPtr<burnup> burnup_(burnup::New(mesh, mat, solverDict));
    
    autoPtr<neutronicsSubSolver> neutronics_
    (neutronicsSubSolver::New(mesh, burnup_, solverDict));
    
    autoPtr<fastFlux> fastFlux_(fastFlux::New(mesh, mat, solverDict));
    
    autoPtr<thermalSubSolver> thermal_
    (thermalSubSolver::New(mesh, mat, solverDict));
    
    autoPtr<mechanicsSubSolver> mechanics_
    (mechanicsSubSolver::New(mesh, mat, rheo_(), solverDict));
    
    autoPtr<fissionGasRelease> fgr_
    (fissionGasRelease::New(mesh, mat, solverDict)); 
    
    autoPtr<gapGasModel> gapGas_
    (gapGasModel::New
        (mesh, mat, thermal_->T(), mechanics_->DorDD(), fgr_,  solverDict)); 

    autoPtr<corrosion> corrosion_
    (corrosion::New(mesh, mat, solverDict)); 

    autoPtr<elementTransport> elementTransport_
    (elementTransport::New(mesh, mat, solverDict));

#ifdef OPENFOAMFOUNDATION
    //- Flow Solver
    autoPtr<flowSubSolver> flow_ = flowSubSolver::New(runTime, solverDict);
#endif
    
    // Check user parameters
    listUserParameters();
    listRegisteredUserParameters();
    checkUserParameters();

    // If post-processing is requested, post-process and then exit
    if (argList::postProcess(argc, argv))
    {
        postProcess(args, runTime);
        return 0;
    } 

    // Keep track of total outer iterations and time steps
    // from start to end
    scalar totalIterations(0);
    scalar totalTimeSteps(0);

    runTime.setEndTime
    (
        min
        (
            heatSrc_->lastTimeMarker(),
            runTime.endTime().value()
        )
    );

    adjustableTime.setLastTimeMarker(runTime.endTime().value());

#ifdef isCommDataLayerIncluded
    const bool isSolveFMI(runTime.controlDict().lookupOrDefault("solveFMI", false));
    fmuControl* fmu = nullptr;
    if (isSolveFMI) fmu = new fmuControl(runTime);
#endif

    // Start time-loop
    while (runTime.run())
    {
#ifdef isCommDataLayerIncluded
        if (runTime.timeIndex() == 0 && isSolveFMI)
        {
            fmu->receive();
            fmu->send();
        }
#endif

        // True when time step is converged
        bool converged(false);

        // Read from fvSolution dict the maximum number of outer iterations
        const dictionary& stressControl = 
        mesh.solutionDict().subDict("stressAnalysis");
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


        runTime++;

        Info<< "Time = " << runTime.userTime() << runTime.unit() << nl << endl;

#ifdef isCommDataLayerIncluded
        do
        {
            if (isSolveFMI) fmu->receive();
#endif

        // Track number of outer iterations (per time-step)
        int nOuterIter = 0;

        // Update mesh if necessary
        if
        (
            mesh.foundObject<fvMesh>("referenceMesh")
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
            Info << "OuterIteration n. " << nOuterIter << nl <<  endl;

            storeGlobalFieldsPrevIter();

            // First update corrosion, only then update AMI. Necessary because
            // if corrosion moves the mesh, the AMI is also updated given that 
            // moving the mesh clears the geometry and the AMIPtr. However this
            // first update, does not use the updated mesh location
            // TODO: find a way to update mesh without AMI?
            corrosion_->correct();

            if(nOuterIter < 1)
            {
                forAll(mesh.boundary(), patchi)
                {             
                    if(isType<regionCoupledOFFBEATFvPatch>(mesh.boundary()[patchi]))
                    {
                    const regionCoupledOFFBEATFvPatch& patch
                        = refCast<const regionCoupledOFFBEATFvPatch>(mesh.boundary()[patchi]);

                        patch.updateAMI();
                    }
                }
            }

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

            nOuterIter++;
            totalIterations++;
            
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
           (nOuterIter < minOuterIter)  
           || 
           (not(converged) && nOuterIter < maxOuterIter) 
        );

#ifdef isCommDataLayerIncluded
        if (isSolveFMI)
        {
            commDataLayer& data = commDataLayer::New(runTime);
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

        // Write fields on disk if this is a write time
        if(adjustableTime.write())
        {
            if (runTime.outputTime())
            {  
                mechanics_->writeStressFields();    
            }
            
            runTime.write();
        }

        //- Check if failure occurred
        mat_->checkFailure();
        
        if (runTime.value() < runTime.endTime().value())
        // TODO Could we have the deltaT setting inside the time handling class?
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
                                << tab << "Current time: " << helperFuns::gather(runTime.value()) << endl
                                << abort(FatalError);
                        }
                    }
                }
            }
        #endif 

            adjustableTime.setDeltaT(deltaT);

#ifdef isCommDataLayerIncluded
            // Last, after all the other setDeltaT
            if (isSolveFMI) fmu->setDeltaT();
#endif
        }

        // Increment time step count
        totalTimeSteps++;

        Info<< "Execution Time = " << runTime.elapsedCpuTime() << " s"
            << "  Clock Time = " << runTime.elapsedClockTime() << " s"
            << "  Total Iterations = " << totalIterations
            << "  Total Timesteps = " << totalTimeSteps
            << nl << endl;
    }

    runTime.writeNow();
    
    clearGlobalFields();
    Info<< "\n end \n";

    return 0;
}


// ************************************************************************* //

