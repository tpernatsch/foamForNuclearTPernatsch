Solver configuration
====================

The thermo-mechanical equations and auxiliary physics solved in a simulation
are defined through the ``solverDict`` dictionary.

This file specifies:

* which physics are active,
* which models are selected and their parameters,
* the mapping between mesh regions and the corresponding material
  definitions, including the property models assigned to each region.

Although materials are declared in this dictionary, their behaviour and
available models are documented separately in :doc:`../materials/index`.

-------------------------------------------------------------------------------

.. rubric:: Location of ``solverDict``

The location of the dictionary depends on how the module is used.

* **Standalone (OFFBEAT)**::

    constant/solverDict

* **Multi-physics (GeN-Foam)**::

    constant/<regionName>/solverDict

Apart from the folder path, the structure and keywords are identical.

-------------------------------------------------------------------------------

.. rubric:: General organization

Each physics is configured through its own subdictionary.
Within that block, the entry ``type`` selects the runtime model, while
additional keywords provide the parameters required by that model.

If a physics dictionary is not present, the corresponding physics is
assumed to be deactivated.

.. note::

   At present, the thermal and mechanical solvers cannot be fully disabled.
   If their dictionaries are omitted, the module falls back to a ``constant``
   mode: fields are read (if present) and made available, but no governing equation is solved.

A simplified example is shown below.

.. code-block:: cpp

    thermalSolver
    {
        type solidConduction;
    }

    mechanicsSolver
    {
        type smallStrain;
        forceSummary on;
    }

    burnup
    {
        type Lassmann;
    }

    fgr
    {
        type SCIANTIX;
    }

    materials
    {
        fuel
        {
            material UO2;
            conductivityModel UO2Matpro;
        }

        cladding
        {
            material Zircaloy;
        }
    }

-------------------------------------------------------------------------------

.. rubric:: Runtime selection principle

For each physics block, the ``type`` entry corresponds to a model
registered in the library.

At runtime, the framework instantiates the appropriate implementation
based on this keyword.

If a model requires additional inputs, they are specified in the same
dictionary.

-------------------------------------------------------------------------------

.. rubric:: Materials selection

The ``materials`` subdictionary assigns a material definition to each
``cellZone`` of the mesh.

While the mapping is declared in ``solverDict``, the description of
available material types and property models is provided in
:doc:`../materials/index`.

Conceptually:

* solver blocks → define **which equations and global models are active**
* materials → define **how matter behaves locally**

-------------------------------------------------------------------------------
 
The following pages describe the available models for each physics.

.. toctree::
   :maxdepth: 1
   :caption: Physics models

   thermal_solver
   mechanics_solver
   burnup
   heat_source
   fast_flux
   fission_gas_release
   gap_gas
   element_transport
   rheology
   slice_mapper
