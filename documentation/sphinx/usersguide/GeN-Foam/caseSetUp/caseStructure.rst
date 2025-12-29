
.. _caseFolderStructureGF:

=====================
Case Folder Structure
=====================

The folder structure of GeN-Foam case is similar to that of typical multi-region solvers shipped with OpenFOAM (e.g., ``chtMultiRegionFoam``). While more traditional OpenFOAM solvers (e.g., ``pimpleFoam`` or ``simpleFoam``) only have the ``0/``, ``constant`` and ``system`` folders, in GeN-Foam, each of these solvers is further subdivided in different sub-folders, one for every region the user is solving for. In the example below, the user is making use of one fluid region, two fluid regions, and one neutronics region. 


.. note ::

    The name of the regions is arbitrary. The ``regionsDict`` dictionary is used to define which physics is solved within each region (see :ref:`Coupling and time stepping <couplingGF>`).
 
.. note ::

    The fields within each of the ``0/`` and ``constant/`` sub-folders depend entirely on the physics that the user decides to solve for in that region. Often, only some fields are mandatory. See the user guide for the ref:`Physics Modules <modules>` for more details.

.. note::

    Regions in a multi-region case can be:

    - **Adjacent**, exchanging data through shared boundaries.
    - **Overlapping**, coupled via mesh-to-mesh projection methods.
    - **Independent**, with no direct interaction.
    
    Any combination of these configurations is supported, enabling flexible multi-physics coupling



.. code-block:: text

    case/
    │
    ├── 0/                              # Initial conditions for all regions
    │   ├── fluidRegion/                # Fluid region fields
    │   │   ├── T                       # Temperature
    │   │   ├── U                       # Velocity
    │   │   └── p_rgh                   # Pressure
    │   ├── solidRegion1/               # Solid region 1 fields
    │   │   ├── T
    │   │   └── D                       # Displacement
    │   ├── solidRegion2/
    │   │   ├── T
    │   └── neutronicsRegion/           # Neutronics region fields
    │       ├── defaultFlux             # Neutron flux
    │
    ├── constant/
    │   ├── fluidRegion/
    │   │   ├── polyMesh/
    │   │   └── g
    │   │   └── phaseProperties
    │   │   └── thermophysicalProperties
    │   │   └── turbulenceProperties
    │   ├── solidRegion1/
    │   │   ├── polyMesh/
    │   │   └── thermoMechanicalProperties
    │   ├── solidRegion2/
    │   │   ├── polyMesh/
    │   │   └── thermoMechanicalProperties
    │   └── neutronicsRegion/
    │       ├── polyMesh/
    │       ├── nuclearData           # Nuclear cross-section data
    │       └── neutronicsProperties    # Solver-specific neutronics settings
    │   └── regionProperties            # Lists all regions 
    │   └── regionsDict                 # Coupling between regions
    │
    ├── system/
    │   ├── fluidRegion/
    │   │   ├── fvSchemes
    │   │   ├── fvSolution
    │   ├── solidRegion1/
    │   │   ├── fvSchemes
    │   │   ├── fvSolution
    │   ├── solidRegion2/
    │   │   ├── fvSchemes
    │   │   ├── fvSolution
    │   └── neutronicsRegion/
    │       ├── fvSchemes
    │       ├── fvSolution
    │   ├── decomposeParDict            # For parallel runs
    │   ├── controlDict                 # Solution control
    │
    ├── postProcessing/
    │
    ├── Allrun
    └── Allclean
``


Folder Overview
===============

At the start of the simulation, a GeN-Foam case folder must contain the following three sub-folders:

1. **``0/`` (or different starting-time folder)**  
   Contains initial and boundary conditions for all fields of all regions.

2. **``constant/``**  
   Contains mesh information (within ``polyMesh/`` ), material properties, behavioral models and correlations,  and module-specific settings.

3. **``system/``**  
   Contains control and numerical settings for the simulation, including discretization schemes (``fvSchemes``), linear solvers (``fvSolution``), and the main simulation controls (``controlDict``). Other files that can be found in this solver include ``blockMeshDict`` (input dictionary for the ``blockMesh`` utility for mesh generation) and ``changeDictionaryDict`` (instructions to modify patch boundaries or other dictionaries via the ``changeDictionary`` utility).

After the simulation is performed, additional folders and files may appear based on the settings configured by the user:

- **Time folders**: Generated during the simulation, they store the field values for specific time steps (e.g., ``0.5/``, ``1.0/``).
- **``postProcessing/`` folder**: Contains data from any post-processing operations performed during or after the simulation, assuming the user has set up at least one ``functionObject``.

.. note::
   The initial folder can also differ from ``0/``. In that case, ensure that either ``fromLatestTime`` is selected as ``startFrom`` in the ``controlDict`` or that ``startTime`` is properly set and ``startFrom`` is configured to ``startTime``.


--------------------------
Optional Files
--------------------------

A GeN-Foam folder may also include optional files, such as:

- **``Allrun`` and ``Allclean`` scripts**: Bash scripts for running and cleaning up cases.
- **``Residuals.gp`` or similar**: Gnuplot scripts for visualizing simulation residuals.
- **``plot.py`` or similar**: Python scripts for plotting data from the ``postProcessing/`` folder.

.. note::
   The ``Allrun`` and ``Allclean`` scripts, as well as plotting scripts, are not directly part of OpenFOAM and are not necessary to successfully run a GeN-Foam simulation.
