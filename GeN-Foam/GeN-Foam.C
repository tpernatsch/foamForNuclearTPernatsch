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

Application
    GeN-Foam

Description
    Multi-physics solver for nuclear reactor analysis. It couples a multi-scale
    fine/coarse mesh 3-phase (liquid, vapour, porous substructure) sub-solver
    for thermal-hydraulics, various sub-solvers for neutronics,
    a displacement-based sub-solver for thermal-mechanics. The
    thermal-hydraulic sub-solver consists of the custom developed FFSEulerFoam
    solver  (https://gitlab.com/virmodoetiae/FFSEulerFoam). It is capable of
    modelling single and two-phase flows, while modelled fuel types consist
    of either liquid fuel (e.g. MSRs) or fuel pin lattices. For the latter,
    the energy dynamics is represented via a 1.5-D finite difference model.

    Reference publications:

    NOTE: these publications do not cover recent multi-phase developments

    Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk, 2015
    "GeN-Foam: a novel OpenFOAM® based multi-physics solver for 2D/3D transient
    analysis of nuclear reactors", Nuclear Engineering and Design 294, pp.
    24-37

    Carlo Fiorina, Konstantin Mikityuk, " Application of the new GeN-Foam
    multi-physics solver to the European Sodium Fast Reactor and verification
    against available codes", Proceedings of ICAPP 2015, May 03-06, 2015 -
    Nice (France), Paper 15226

    Authors of this file (and associated .C or included .H files):
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>
    Stefan Radman <stefanradman92@gmail.com; stefan.radman@epfl.ch;>
    EPFL (Switzerland)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "regionSolvers.H"
#include "setDeltaT.H"
#include "meshToMesh.H"

// #include <gperftools/profiler.h>

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

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // ProfilerStart("gen-foam.prof");  // Start profiling

    argList::addBoolOption
    (
        "legacy",
        "Run GeN-Foam in legacy mode (V.1)"
    );

    argList::addBoolOption
    (
        "initializeMappedFields",
        "Map fields across regions as specified in multiRegionCouplingDict"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const bool legacy(args.found("legacy"));
    const bool mappingMode(args.found("initializeMappedFields"));

    if (legacy)
    {
        Info<< nl << "Running GeN-Foam in legacy mode (V.1)" << nl << endl;
        #include "modifyDictLegacy.H"

        Info<< nl << "Dictionaries modified" << nl << endl;
    }

    // Create the region meshes and solvers
    regionSolvers solvers(runTime);
    forAll(solvers, i)
    {
        solvers[i].correctBaffleLessFields();
    }

    // Set the initial time-step
    setDeltaT(runTime, solvers);


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Starting time loop\n" << endl;

    solvers.setGlobalPrefix();

    solvers.mapper().initializeMappedFields(runTime, mappingMode);

    #ifdef isCommDataLayerIncluded
    // FMU solution control
    const bool isSolveFMI(runTime.controlDict().lookupOrDefault("solveFMI", false));
    fmuControl* fmu = nullptr;
    if (isSolveFMI) fmu = new fmuControl(runTime);
    #endif

    while (runTime.run() and !(mappingMode))
    {
        #ifdef isCommDataLayerIncluded
        // At the first time step
        if (runTime.timeIndex() == 0 && isSolveFMI)
        {
            fmu->receive();
            fmu->send();
        }
        #endif

        solvers.setGlobalPrefix();

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        solvers.mapper().mapAllFields(runTime);

        #ifdef isCommDataLayerIncluded
        do
        {
        if (isSolveFMI) fmu->receive();
        #endif

        // Solve each physics once (loose coupling)

        forAll(solvers, i)
        {
            solvers[i].deformMesh();
            solvers[i].correctPhysics();
            solvers[i].correctBaffleLessFields();
        }

        #ifdef isCommDataLayerIncluded
        if (isSolveFMI)
        {
            // Update isConverged flag as FMU output
            commDataLayer& data = commDataLayer::New(runTime);
            scalar& isConvergedFMI = data.getObj<scalar>
            (
                "isConverged",
                commDataLayer::causality::out
            );
            isConvergedFMI = 1.0;

            // forAll(solvers, i)
            // {
            //     if (solvers[i].getResidual() > 1e-4)
            //     {
            //         isConvergedFMI = 0.0;
            //     }
            // }
        }

        if (isSolveFMI) fmu->send();
        } // End fmu implicit loop
        while (isSolveFMI && fmu->loop());
        #endif 

        // Adjust the time-step according to the solver maxDeltaT
        // Need to be after the FMI loop
        adjustDeltaT(runTime, solvers);

        #ifdef isCommDataLayerIncluded 
        // Last, after all the other setDeltaT
        if (isSolveFMI)
        {
            fmu->setDeltaT();
            if (fmu->isTerminated()) break;
        }
        #endif

        solvers.setGlobalPrefix();

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    // ProfilerStop();  // Stop profiling

    return 0;
}


// ************************************************************************* //
