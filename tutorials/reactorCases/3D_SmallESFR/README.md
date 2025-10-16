# 3D Small ESFR

3D_SmallESFR is a slightly simplified 3-D model of the core of the European Sodium Fast Reactor (ESFR). It is a full steady-state multiphysics case, including thermal-hydraulics, eigenvalue neutronics and thermal-mechanics (with expansion of the meshes based on the displacement field calculated by the thermal-mechanics). 
This tutorial also include an example of Gmsh scripts for generating reactor geometries. The scripts can be found in the folder 3dESFRMesh. In particular, 4 folders are present, each including a gmsh script for: neutronics; thermal-mechanics; thermal-hydraulics; thermal-hydraulics with one separate boundary conditions for one of the assemblies, allowing to simulate the case of a blocked assembly. Meshes can be created by running the command line "gmsh esfr.main". This will open gmsh (that must be first installed on your computer) and, after quite a long time, show you a geometry. This geometry must be meshed and the resulting mesh saved as a .msh file. These .msh file can be converted to an OpenFOAM-readable format (the polymesh folder) using the gmshToFoam tool distributed with OpenFOAM. This should be done at the root of the case folder and it will create a polyMesh folder in the "constant" folder. Depending on the mesh that is being created, this folder must then be copied in the "constant" folder of the corresponding physics. Please notice that the mesh generation is not a necessary step for running the tutorials. The polyMesh folders are already present.

Two cases are present, to show the usage of two different thermomechanical solvers, i.e. the legacy linear-elastic solver and the newly implemented advanced thermomechanics solver

The steady state solutions of the legacy approach and the new one are compared. The `./Allrun_parallel` script runs both cases sequentially on 4 cores.

The same geometry is employed in a guided tutorial with partially filled dictionaries in [ESFRWorkshop](./../../guidedCases/ESFRWorkshop/).
