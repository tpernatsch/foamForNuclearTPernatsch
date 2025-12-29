.. _meshingGF:

Meshing
=======

GeN-Foam uses any number of meshes for that can be e.g neutronics,
thermal-hydraulics and thermal-mechanics. There is no requirement for the meshes
to occupy the same region of space. Consistent mapping of fields is performed
and a reference value is given to a field if no correspondence is found in the
mesh where its value is being projected from. It follows that the geometry for
neutronics can cover only a small part of the overall reactor geometry. Meshes
can be created with every OpenFOAM-compatible tool (e.g blockMesh, ANSYS/FLUENT,
...). Meshes can (should) be divided into zones (cellZones) to allow the use of
different physical properties (e.g., cross-sections, power models, ...) in
different reactor regions. Sometimes, when converting a mesh to the OpenFOAM
(``polyMesh`` folder) format, cellSets (and not cellZones) are created. The
``topoSetDict`` can be used to convert cellSets to cellZones.

The `3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/>`_
tutorial includes an example of mesh generation with Gmsh `GMSH2009 <https://onlinelibrary.wiley.com/doi/10.1002/nme.2579>`_.
The folder includes the subfolder *3dESFRMesh* containing three subfolders for the
generation of the meshes for neutronics, thermal-hydraulics and
thermal-mechanics. To create the 3D ESFR core geometry, one should execute the
following command in a terminal:

.. code :: bash

    gmsh esfMain.geo


This will create (after a relatively long time) a geometry and open the Gmsh
graphical interface. The geometry must then be meshed and the resulting mesh
saved and copied in the root of the *3D_SmallESFR* folder.

By typing in a terminal:

.. code :: bash

    gmshToFoam nameMeshFile


A *polyMesh* folder will be created (or updated) in the folder *constant*. One
should then copy this folder in *constant/neutroRegion*, *constant/fluidRegion*
or *constant/thermoMechanicalRegion*, and repeat the operation for all meshes.

Please notice that the *3D_SmallESFR* tutorial already contains the correct
*polyMesh* folders so that one can avoid the mesh generation step.