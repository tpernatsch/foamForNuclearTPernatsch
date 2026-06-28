# 3D Small ESFR

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-onePhase-blue.svg)]() [![badge](https://img.shields.io/badge/Neutronics-diffusion-blue.svg)]() [![badge](https://img.shields.io/badge/ThermalMechanics-extendedThermoMechanics-blue.svg)]() [![badge](https://img.shields.io/badge/Multiphysics-looseCoupling-orange.svg)]()


## Description

`3D_SmallESFR` is a slightly simplified 3-D model of the core of the European Sodium Fast Reactor (ESFR). It is a full steady-state multiphysics case, including thermal-hydraulics, eigenvalue neutronics, and thermal-mechanics (with mesh expansion based on the displacement field computed by the thermal-mechanics).

This tutorial also includes an example of Gmsh scripts for generating reactor geometries. The scripts are located in the `3dESFRMesh` folder. In particular, four folders are provided, each containing a Gmsh script for: neutronics, thermal-mechanics, thermal-hydraulics, and thermal-hydraulics with one separate boundary condition for one of the assemblies. This last option allows you to simulate the case of a blocked assembly.

Meshes can be created by running the command line `gmsh esfr.main`. This will open Gmsh (which must be installed on your computer) and, after a fairly long time, display the geometry. The geometry must then be meshed, and the resulting mesh saved as a `.msh` file.

These `.msh` files can be converted to an OpenFOAM-readable format (the `polyMesh` folder) using the `gmshToFoam` tool distributed with OpenFOAM. This should be done at the root of the case folder, and it will create a `polyMesh` folder in the `constant` folder. Depending on the mesh being created, this `polyMesh` folder must then be copied into the `constant` folder of the corresponding physics.

Please note that mesh generation is not required to run the tutorials. The `polyMesh` folders are already provided.

Two cases are included to demonstrate the use of two different thermomechanical solvers: the legacy linear-elastic solver and the newly implemented advanced thermomechanics solver.

The steady-state solutions obtained with the legacy approach and the new one are compared. The `./Allrun_parallel` script runs both cases sequentially on 4 cores.

The same geometry is also used in a guided tutorial with partially filled dictionaries in [ESFRWorkshop](./../../guidedCases/ESFRWorkshop/).