.. _meshingGF:


=========================
Meshing in foamForNuclear
=========================

Meshing is a fundamental step in setting up any OpenFOAM-based simulation, including multi-physics cases in foamForNuclear. This section introduces general concepts, region-specific requirements, and workflows for external meshing tools.


General Concepts
================

OpenFOAM uses **polyhedral meshes**, which are stored in the ``constant/polyMesh`` folder of each region. A mesh consists of points, faces, and cells, organized in a way that supports finite-volume discretization.

Meshes can be generated using:
- **OpenFOAM utilities** such as ``blockMesh`` (structured meshes) or ``snappyHexMesh`` (unstructured meshes).
- **External tools** like Salome, Gmsh, Cubit, various ANSYS tools, etc.
- **Conversion utilities**: OpenFOAM ships with many mesh converters (e.g., ``fluentMeshToFoam``, ``gmshToFoam``, ``cubitToFoam``) to import meshes from other platforms.

.. note::
    Each region in a multi-region case (e.g., fluid, solid, neutronics) requires its own mesh stored under ``constant/<regionName>/polyMesh``.


Cell Zones
==========

Many physics modules in foamForNuclear require **cellZones** to define subsets of the mesh for applying models (e.g., fuel regions, cladding, coolant channels). A **cellZone** is a named collection of cells within a region.

- In **Salome**, cellZones correspond to **groups of volumes** created during geometry partitioning.
- In **Gmsh**, they are defined as **physical volume groups**.
- In **Cubit**, zones are typically created as **block sets**.

These groups are preserved during mesh export and recognized by OpenFOAM converters, which create entries in ``constant/<regionName>/cellZones``.


External Meshing Tools
======================

foamForNuclear supports meshes from both proprietary and open-source tools. Below are exemplary workflows:

**1. Salome Workflow**
- Create geometry and partition it into volumes.
- Assign **groups** for each physical region (fuel, cladding, coolant).
- Generate the mesh and export in UNV format.
- Convert to OpenFOAM using ``ideasUnvToFoam``.
- Verify that cellZones are correctly imported.

**2. Gmsh Workflow**
- Build geometry and define **Physical Groups** for volumes.
- Generate the mesh and export in Gmsh format (.msh).
- Convert using ``gmshToFoam``.
- Check ``constant/<regionName>/cellZones`` for zone definitions.

.. tip::
    Always check mesh quality using ``checkMesh`` after conversion. Poor-quality cells can lead to solver instability.

.. warning::
    Sometimes, when converting a mesh to the OpenFOAM (``polyMesh`` folder) format, cellSets (and not cellZones) are created. The ``topoSetDict`` utility can be used to convert cellSets to cellZones.


Changing Patch Type (Optional)
==============================

Boundary types (e.g., ``wedge``, ``empty``, ``coupled``) may need adjustment after mesh generation. Use:

- **``foamDictionary``**: Reads ``system/changeDictionaryDict`` to modify patch definitions or other entries.
  
.. warning::
    ``foamDictionary`` has replaced ``changeDictionary`` that was used in older OpenFOAM versions.




Example
=======

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


Summary
=======

- Each region requires its own mesh in ``constant/<regionName>/polyMesh``.
- CellZones are essential for applying material properties and physics models.
- OpenFOAM provides robust utilities for mesh generation and conversion.
- Salome, Gmsh, Cubit and Fluent are widely used for generating complex geometries in OpenFOAM.