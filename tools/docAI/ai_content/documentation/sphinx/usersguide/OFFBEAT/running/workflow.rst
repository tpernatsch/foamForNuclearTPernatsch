Workflow
========

OFFBEAT operates similarly to typical solvers shipped with OpenFOAM (e.g. ``icoFoam``,
``pisoFoam``). The user provides a mesh, a control dictionary (``controlDict``), solution
parameters (``fvSolution``), discretization schemes (``fvSchemes``), an OFFBEAT-specific
solver dictionary (``solverDict``), and initial and boundary conditions for the main fields
in the initial time-step folder (e.g. ``0/``).

This guide provides an overview of the main steps and commands required to run OFFBEAT.
Each part of the workflow is described in more detail in the following sections of this
User’s Guide. When executed sequentially from the case folder, these steps define a
complete and consistent procedure for running the code.


List of main steps and commands
-------------------------------


Building the mesh
-----------------

OpenFOAM provides several options for mesh generation. Meshes can be created using external
tools and converted to OpenFOAM format, or generated directly within OpenFOAM using dedicated
utilities. For one- and two-dimensional r–z simulations (though not necessarily for
r–\ :math:`\theta`), the meshes typically used for nuclear fuel simulations are relatively
simple block-structured meshes. In this case, OpenFOAM’s built-in ``blockMesh`` utility is
often sufficient.

``blockMesh`` builds the computational mesh based on the ``blockMeshDict`` dictionary located in the
``system/`` directory. Executing this command creates the ``constant/polyMesh/`` folder, which stores all mesh information
(points, faces, cells, etc.) in a set of dictionaries.

  .. note::

    OFFBEAT provides custom tools to automatically generate ``blockMeshDict`` files for
    standard fuel-rod geometries. Further details are provided in the Tools section of this
    User’s Guide.


Changing patch types (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Whether the mesh is generated using ``blockMesh`` or imported from external software, the
boundary patch types (e.g. wedge, empty, coupled) may not always be correctly defined.
Similarly, a mesh provided by a third party may require adjustments to patch or field
definitions. In such cases, the ``changeDictionary`` utility can be used to modify dictionaries
and files in the case folder. By default, it reads instructions from ``system/changeDictionaryDict``.
In OFFBEAT test cases, this utility is sometimes used to modify entries in the ``polyMesh/boundary``
dictionary, or to adjust patch definitions and field settings.

  .. warning::

    Although still used in OFFBEAT test cases, the ``changeDictionary`` utility is
    deprecated in recent OpenFOAM versions and has been superseded by ``foamDictionary``.


Setting up the simulation
-------------------------

Once the mesh is ready, the simulation can be configured by defining the following
components:

1. **Boundary and initial conditions**

   Defined in the ``0/`` folder (or another starting-time folder). These typically include
   temperature (``T``), displacement (``D``), and gas gap composition (``gapGas``).

2. ``solverDict``

   Central dictionary specifying the physical models, numerical options, and material
   properties.

3. ``system/`` **settings**

   Configuration of ``controlDict`` (start and end times, writing intervals, function
   objects), ``fvSchemes`` (discretization schemes), and ``fvSolution`` (solver types and
   convergence tolerances).


Running the code
----------------

The command ``offbeat`` in the terminal executes OFFBEAT. Run-time information (the log) is printed
to the terminal. In most cases, the successful completion of the simulation is indicated by the
presence of the ``end`` string at the end of the output. This should always be checked.

Several options are commonly used to handle the log output:

  * ``offbeat | tee log.offbeat``
    Writes the output both to the terminal and to a log file. The name ``log.offbeat`` is
    recommended for compatibility with auxiliary scripts and utilities.

  * ``offbeat | tee log.offbeat & gnuplot Residuals.gp``
    Runs OFFBEAT while simultaneously launching the ``Residuals.gp`` gnuplot script (if
    present in the case folder) to monitor residual convergence during the simulation.

  * ``offbeat > log.offbeat``
    Redirects the output exclusively to a log file without displaying it in the terminal.


Post-processing (optional)
--------------------------

OpenFOAM provides extensive post-processing capabilities. Most users rely on either
ParaView or OpenFOAM’s built-in post-processing utilities. Additional information is
provided in the Data Processing section of this User’s Guide.

* ``paraFoam``
  Enables visualization of the fields computed by OFFBEAT using ParaView. Only fields
  written at stored time steps can be visualized.

* ``postProcess -func "sampleDict"``
  Example usage of the ``postProcess`` utility to sample fields at selected locations in
  the computational domain. The utility reads instructions from ``system/sampleDict``, and
  results are written to the ``postProcessing/`` folder.

  .. note::

    Only fields written at stored time steps can be sampled or probed in this way. Similar
    functionality can be configured in the ``functions`` subdictionary of ``controlDict``
    to extract data during the simulation.

* ``python plot.py``
  Python scripts can be used to automatically plot quantities such as centreline
  temperature or local burnup from data stored in ``postProcessing/``. Test cases include
  examples that plot the centreline temperature at the rod centre together with the
  corresponding local burnup.


Case cleaning
-------------

* ``foamListTimes -rm``
  Deletes all stored time-step folders except ``0/``, effectively removing the results of
  a previous simulation.

  If available, the ``Allclean`` script can be used instead to remove not only time folders
  but also auxiliary files such as logs and generated figures.