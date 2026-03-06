Thermo-mechanics
========================

This section describes how to configure thermo-mechanical simulations using the
``thermoMechanics`` module of foamForNuclear.

While the module solves the classical conservation of energy and momentum in solids,
it is specifically designed for **nuclear materials applications** and therefore provides
dedicated capabilities that go well beyond standard thermo-mechanics.  
These include, for example:

* burnup evolution,
* fission gas production and release,
* gaseous and solid swelling,
* irradiation-dependent material behaviour,
* gap heat transfer and contact mechanics,
* element and isotope transport.

The information provided here applies both when the module is:

* used as a standalone application (e.g. OFFBEAT), and  
* coupled within a GeN-Foam multi-physics simulation.

The physics models, boundary conditions, numerical controls, and discretization
settings are identical in both situations.

The only practical difference is the case folder structure:

* **Standalone (OFFBEAT)** → files are located directly in ``0/``, ``constant/`` and ``system/``.
* **Multi-physics (GeN-Foam)** → the same files are located under  
  ``0/<regionName>/``, ``constant/<regionName>/`` and ``system/<regionName>/``.

Apart from this change in path, dictionaries and keywords are unchanged.

.. warning::

   This guide assumes familiarity with OpenFOAM case organization and
   dictionary syntax. Users should already be comfortable with mesh
   preparation, boundary condition assignment, and solver execution.

-------------------------------------------------------------------------------

.. rubric:: What is configured in this module?

A thermo-mechanical case requires the definition of:

* the physics and material models,
* the boundary and initial conditions,
* the numerical solution strategy,
* the time management,
* the discretization schemes.

Model and materials selection is primarily handled through the ``solverDict`` file,
while fields and numerical behaviour are defined in ``0/`` and ``system/``.

-------------------------------------------------------------------------------

.. toctree::
   :numbered:
   :maxdepth: 1
   :caption: Contents

   boundary_conditions/index
   solver_configuration/index
   materials/index
   numerics/solution_control
   numerics/schemes