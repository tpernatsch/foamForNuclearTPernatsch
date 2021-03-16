/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    GeN-Foam

Description
    Multi-physics solver for nuclear reactor analysis. It couples a multi-scale
    fine/coarse mesh 3-phase (liquid, vapour, porous substructure) sub-solver 
    for thermal-hydraulics, a multi-group diffusion sub-solver for neutronics,
    a displacement-based sub-solver for thermal-mechanics. The 
    thermal-hydraulic sub-solver consists of the custom developed FFSEulerFoam 
    solver  (https://gitlab.com/virmodoetiae/FFSEulerFoam). It is capable of
    modelling single and two-phase flows, while modelled fuel types consist
    of either liquid fuel (e.g. MSRs) or fuel pin lattices. For the latter,
    the energy dynamics is represented via a 1.5-D finite difference model.

Reference publications
    
    NOTE: these publications do not cover recent multi-phase developments

    Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk, 2015
    "GeN-Foam: a novel OpenFOAM® based multi-physics solver for 2D/3D transient
    analysis of nuclear reactors", Nuclear Engineering and Design 294, pp. 
    24-37

    Carlo Fiorina, Konstantin Mikityuk, " Application of the new GeN-Foam 
    multi-physics solver to the European Sodium Fast Reactor and verification 
    against available codes", Proceedings of ICAPP 2015, May 03-06, 2015 - 
    Nice (France), Paper 15226

Authors of this file (and associated .C or included .H files)
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>
    Stefan Radman <stefanradman92@gmail.com; stefan.radman@epfl.ch;>
    EPFL (Switzerland)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "SquareMatrix.H"
#include "fvMatrixExt.H"
#include "meshToMesh.H"
#include "regionProperties.H"
#include "mergeOrSplitBaffles.H"
#include "volPointInterpolation.H"
#include "fixedGradientFvPatchFields.H"
#include "UPstream.H"

#include "multiphysicsControl.H"
#include "thermalHydraulicModel.H"
#include "neutronics.H"
#include "thermoMechanics.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"
    #include "createMeshInterpolators.H"
    #include "createCouplingFields.H"
    #include "createOutput.H"

    Info<< "\nStarting time loop\n" << endl;

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" 
        << nl << endl;

    #include "setDeltaT.H"

    while (runTime.run())
    {
        runTime++;

        #include "setDeltaT.H"
        
        Info << "Time = " << runTime.timeName() << nl << endl;

        while (multiphysics.loop())
        {
            #include "solve.H"
        }

        runTime.write();

        #include "writeOutput.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
