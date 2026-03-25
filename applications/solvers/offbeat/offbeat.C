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
#include "userParameters.H"
#include "dynamicFvMesh.H"

#include "timeHandling.H"
#include "offbeatTime.H"
#include "extendedThermoMechanics.H"

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

    // Time handling class, mainly for time stepping
    timeHandling adjustableTime(mesh, runTime);

    // Create thermomechanics object
    Foam::solvers::extendedThermoMechanics thermoMech(mesh);
    
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

    // Keep track of total time steps from start to end
    scalar totalTimeSteps(0);

    // Use smallest between controlDict.endTime and last heat source list time point
    runTime.setEndTime(
        min
        (
            thermoMech.endTime(), 
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
        runTime++;

        Info<< "Time = " << runTime.userTime() << runTime.unit() << nl << endl;

        #ifdef isCommDataLayerIncluded
        thermoMech.correctPhysics(fmu);
        #else
        thermoMech.correctPhysics();
        #endif

        // Write fields on disk if this is a write time
        if(adjustableTime.write())
        {          
            runTime.write();
        }
        
        if (runTime.value() < runTime.endTime().value())
        // TODO Could we have the deltaT setting inside the time handling class?
        {   
            scalar deltaT = thermoMech.maxDeltaT();

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
            << "  Total Iterations = " << thermoMech.totalIterations()
            << "  Total Timesteps = " << totalTimeSteps
            << nl << endl;
    }

    runTime.writeNow();
    
    clearGlobalFields();
    Info<< "\n end \n";

    return 0;
}


// ************************************************************************* //

