Case Folder Structure
=====================

The folder structure of an OFFBEAT case is similar to that of typical solvers shipped with
OpenFOAM (for example, ``icoFoam`` and ``pisoFoam``).

.. code-block:: text

    case/
    │
    ├── 0/                     # Initial conditions folder (or starting time folder)
    │   ├── T                  # Temperature field
    │   ├── D                  # Total displacement field
    │   └── gapGas             # Gas gap composition field
    │
    ├── constant/              # Contains mesh and solver-specific settings
    │   ├── polyMesh/          # Mesh information (points, faces, etc.)
    │   ├── solverDict         # Central dictionary specifying models and material properties
    │   └── input_settings.txt # Input file for SCIANTIX (if the fgrSCIANTIX model is activated)
    │
    ├── system/                # Contains control and numerical settings
    │   ├── controlDict        # Simulation control parameters
    │   ├── fvSchemes          # Discretization schemes
    │   └── fvSolution         # Linear solver settings
    │
    ├── postProcessing/        # Generated during the simulation (probes, sampling, etc.)
    │
    ├── Allrun                 # Optional script to run the case
    ├── Allclean               # Optional script to clean the case
    └── input_settings.txt     # SCIANTIX input file (legacy option)


Folder overview
---------------

At the start of a simulation, an OFFBEAT case folder typically contains the following three
subfolders:

1. ``0/`` (or an alternative starting-time folder), which defines the initial and boundary
   conditions for the main solution fields, such as temperature (``T``), displacement
   (``D``), and gas gap composition (``gapGas``).

2. ``constant/``, which contains information that remains unchanged during the simulation,
   including the computational mesh (stored in the ``polyMesh`` subfolder), material
   properties, and solver-specific settings defined in the ``solverDict`` file.

3. ``system/``, which contains the control and numerical settings governing the simulation,
   including ``controlDict``, ``fvSchemes``, and ``fvSolution``.

.. figure:: ../images/case_folder.svg
   :width: 600
   :align: center

   Typical folder structure for an OFFBEAT simulation

After the simulation has been executed, additional folders and files may appear depending on
the user configuration, in particular:

* Additional time folders (e.g. ``0.5/`` and ``1.0/``), generated according to the selected
  write-time options and containing the field values stored at specific time steps.

* The ``postProcessing/`` folder, which contains data produced by OpenFOAM function objects,
  such as probes and sampling utilities.



``0/`` (or different starting-time folder)
------------------------------------------

This folder contains the initialization of the main solution fields.

.. note::

   The initial folder does not necessarily have to be ``0/``. In that case, the user must ensure
   that either ``fromLatestTime`` is selected as ``startFrom`` in ``controlDict``, or that
   ``startTime`` is properly set and ``startFrom`` is configured to ``startTime``.

Typically, the initial folder includes one file per field requiring initial and boundary
conditions. The main fields are:

* Temperature (``T``)
    Primary thermal field solved by the heat conduction model.

* Total displacement (``D``)
    Displacement field used by total Lagrangian mechanical solvers. For incremental
    formulations, the corresponding incremental displacement field (``DD``) is used instead.

* Gas gap composition (``gapGas``)
    Dictionary defining the initial gas composition and pressure in the fuel–cladding
    gap.

.. warning::

   For backward compatibility, ``gapGas`` may be located directly in the ``0/`` folder.
   However, it is strongly recommended to place it in ``0/uniform/`` so that it is correctly
   distributed in parallel simulations. After each time step, OFFBEAT writes the updated
   ``gapGas`` dictionary to ``time_step_number/uniform/``.

.. figure:: ../images/0_folder.svg
   :width: 600
   :align: center

   Contents of the ``0/`` folder. The recommended location for ``gapGas`` is ``0/uniform/``.


``constant/`` folder
--------------------

This folder contains information that remains constant throughout the simulation.

Key elements include:

* The ``polyMesh/`` folder, which stores the mesh connectivity information (points, faces,
  owner/neighbour lists). It is typically created using mesh generation utilities such as
  ``blockMesh``.

* The ``solverDict`` file, which serves as the central input dictionary for OFFBEAT and
  defines the physical models, numerical options, and material properties used in the
  simulation.

  Most OFFBEAT-specific input is concentrated in ``solverDict``. In contrast to standard
  OpenFOAM solvers, where configuration is distributed across multiple dictionaries, OFFBEAT
  adopts a more centralized and consolidated input structure.

.. figure:: ../images/constant_folder.svg
   :width: 600
   :align: center

   Contents of the ``constant/`` folder


``system/`` folder
------------------

This folder contains the main control dictionaries governing the execution of the simulation.

Typical contents include:

* ``blockMeshDict``, the input dictionary for the ``blockMesh`` utility, defining vertices,
  blocks, and boundary patches for structured mesh generation.

    .. note::

        ``blockMeshDict`` is required only when the ``blockMesh`` utility is used. Alternatively,
        meshes can be generated using external tools such as Salome or Gmsh and imported into
        OpenFOAM using the available mesh conversion utilities.

* ``changeDictionaryDict`` (optional), the input dictionary used by the ``changeDictionary``
  utility to modify patch definitions or other dictionaries in the case setup.

    .. warning::

        Although still used in several OFFBEAT test cases, the ``changeDictionary`` utility is
        deprecated in recent OpenFOAM versions and has been superseded by ``foamDictionary``.

* ``controlDict``, which defines the main run-time controls, including start and end times,
  time-step management, write intervals, and configured function objects.

* ``fvSchemes``, which specifies the discretization schemes used for time derivatives,
  gradients, divergences, and other differential operators.

* ``fvSolution``, which defines the numerical solver settings, including linear solver selection,
  tolerances, and convergence criteria.


Optional files
--------------

An OFFBEAT case may also include optional auxiliary files, depending on the workflow and
post-processing needs:

* ``input_settings.txt``
    Input file for the SCIANTIX fission gas release module. In recent OFFBEAT versions, SCIANTIX
    settings are preferably embedded directly in the ``fgrOptions`` subdictionary of
    ``constant/solverDict``. For backward compatibility, if these settings are not found there,
    OFFBEAT will look for this file in the root case directory.

* ``Allrun`` and ``Allclean``
    Shell scripts provided for convenience to automate case execution and cleanup.

* ``Residuals.gp`` or similar
    Gnuplot scripts used to visualize residual histories during or after a simulation.

* ``plot.py`` or similar
    Python scripts for custom post-processing of data stored in the ``postProcessing/`` folder.

.. note::

    The ``Allrun`` and ``Allclean`` scripts, as well as plotting utilities, are not required by
    OpenFOAM or OFFBEAT to successfully run a simulation. They are provided purely for user
    convenience.