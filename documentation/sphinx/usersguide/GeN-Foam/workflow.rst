
.. _workflowGF:

==========================
Workflow for GeN-Foam
==========================

GeN-Foam is a multi-region, multiphysics solver built on OpenFOAM. It operates similarly to multi-region solvers of OpenFOAM. Each region (e.g., fluid, solid, neutronics) requires its own mesh and dictionaries, organized under the top-level folders ``0/``, ``constant/``, and ``system/``.

This section provides an overview of the main steps and commands for setting up and running GeN-Foam. Each part of the workflow is detailed in subsequent sections of this User Guide.

Designing Your Simulation Case
==============================


Before running GeN-Foam, the user must carefully design the simulation case. This involves defining the physical problem and translating it into a multi-region setup. Key decisions include:

- **Number of regions**: Identify which domains are needed (e.g., fluid, solid, neutronics).
- **Physics modules per region**: Identify the appropriate :ref:`physics module <modules>`  for each region (e.g., thermal-hydraulics for fluid, thermal-mechanics for structural components, neutron transport for the core) and, when needed, the specific model within the module (e.g., single-phase thermal-hydraulics, linear eleasticity, and diffusion).
- **Coupling strategy**: Determine how regions interact—through shared boundaries, mesh-to-mesh projections, or act indepndently.
- **Material and operational data**: Prepare thermophysical properties, nuclear cross-sections, etc.

A well-structured case ensures accurate multiphysics coupling and efficient simulation. Poor planning can lead to inconsistencies in mesh interfaces, missing dictionaries, or incorrect solver configurations.


.. note::

    In GeN-Foam, one region is always associated with a single module. For multi-physics simulations involving multiple modules on the same region of space (e.g., porous-medium thermal-hydraulics and diffusion neutronics in a reactor core), the user must employ multiple overlapped meshes coupled via mesh-to-mesh projection. When the two meshes are identical, this approach provides the same results as a single-mesh approach, with minimal computational overhead.  



List of Main Steps and Commands
===============================


Building the Mesh(es)
---------------------

OpenFOAM uses **polyhedral meshes**, stored in ``constant/<regionName>/polyMesh`` for each region. Meshes can be generated using:

- **OpenFOAM utilities**: ``blockMesh`` for structured meshes, ``snappyHexMesh`` for complex geometries.
- **External tools**: e.g., Salome, Gmsh, Cubit.
- **Conversion utilities**: OpenFOAM provides converters such as ``gmshToFoam``, ``ideasUnvToFoam``, and ``cubitToFoam``.

.. note::
    Each region (fluid, solid, neutronics) must have its own mesh. GeN-Foam requires consistent interfaces between regions for coupling.



Setting Up the Simulation
-------------------------

After meshing, configure:

- **``constant/regionProperties``**: Lists all regions.
- **``constant/regionsDict``**: Set the coupling among regions.
- **Region-specific dictionaries**:
  - ``constant/<regionName>/dictionaryName``: set properties, behavioral models, etc., for each region.
- **Initial and boundary conditions**: Defined in ``0/<regionName>/`` for fields like:
  - Fluid: ``T``, ``U``, ``p``.
  - Solid: ``T``, ``D``.
  - Neutronics: ``flux``.
- **System settings**:
  - ``system/<regionName>/controlDict`` for time control.
  - ``fvSchemes`` and ``fvSolution`` for discretization and solvers.


Running the Code
----------------

- **``GeN-Foam``**: Executes the multiphysics solver. Monitor output in the terminal or redirect to a log file:

  - ``GeN-Foam | tee log.genFoam``: Output to terminal and log file.
  - ``GeN-Foam > log.genFoam``: Output only to log file.
  - ``mpirun -np 8 GeN-Foam -parallel | tee log.genFoam``: Run in parallel on 8 cores. Output to terminal and log file.

.. tip::
    Check for the ``End`` string in the log to confirm successful completion.

.. warning::
    Do not forget``-parallel`` when running in parallel!

Post-Processing
---------------

Visualization and data extraction options:

- **``paraFoam``**: Launches ParaView for visualization of written time steps (requires ParaView).
- **``touch para.foam``** followed by **``paraview para.foam``**: Launches ParaView for visualization of written time steps in case paraFoam was not installed (requires ParaView).
- **``postProcess -func sampleDict``**: Samples fields based on ``system/sampleDict``; results stored in ``postProcessing/``.
- **Python scripts**: Automate plotting of quantities like temperature profiles or power distribution.

.. note::
    Only fields from written time steps can be visualized or sampled. Finer results can be stored in ``postProcessing/`` by using function objects.


Case Cleaning
-------------

- **``foamListTimes -rm``**: Removes all time-step folders except ``0/``.
- **``Allclean``**: This bash script typically provided with eacch tutorial. It deletes time folders and auxiliary files (logs, PNGs, etc.).


