Materials
=============================

In OFFBEAT, all material-related information is defined in the
``materials`` dictionary of the ``solverDict``.

Materials are assigned **by mesh region**, using ``cellZone`` objects
defined during the meshing stage. Each material entry specifies:

* the **material type** (e.g. UO₂, MOX, Zircaloy),
* the **material properties** (e.g. thermal conductivity, density),
* optional **behavioral models** (e.g. swelling, densification),
* **constitutive behavior** (e.g. elasticity, plasticity),
* optional **damage models**.

-------------------------------------------------------------------------------

.. rubric:: What is a ``cellZone``

In OpenFOAM, a ``cellZone`` is a named group of mesh cells.
It is used to identify regions of the domain that share common physical
properties or modeling assumptions.

``cellZone`` objects are:

* defined **at meshing time** (e.g. via ``blockMesh`` or ``topoSet``),
* stored in the mesh and remain fixed throughout the simulation,
* commonly used to distinguish different materials or structural parts
  (e.g. fuel, cladding, moderator).

-------------------------------------------------------------------------------

.. rubric:: Material assignment by ``cellZone``

The ``materials`` dictionary must contain **one subdictionary per
``cellZone``** defined in the mesh.

The following rules apply:

* The name of each subdictionary **must exactly match** a ``cellZone`` name.
* Each ``cellZone`` must appear **once and only once**.
* No ``cellZone`` may be omitted.
* Multiple ``cellZone`` objects may use the **same material type**, but
  each must still have its own entry.
* If several ``cellZone`` objects share the **same material type and
  options**, they may be grouped using **wildcards**, provided the
  resulting pattern uniquely matches all intended zones.

The ``cellZone`` name is purely geometric and has **no physical meaning**
by itself within the simulation.

-------------------------------------------------------------------------------

.. rubric:: Material type selection

The physical nature of the material (e.g. UO₂, MOX, Zircaloy) is defined
by the ``type`` keyword inside each material subdictionary.

.. code-block:: cpp

   materials
   {
       fuel
       {
           type UO2;
       }
   }

* The ``type`` keyword selects the material class.
* It determines which properties, behavioral models, and constitutive
  laws are available for that material.

.. warning::

   The deprecated keyword ``material`` is still accepted for backward
   compatibility but should no longer be used.

-------------------------------------------------------------------------------

.. rubric:: Property model selection

Material properties are configured through dedicated property
subdictionaries. Each property defines the model used to compute it via
the ``type`` keyword.

.. code-block:: cpp

   materials
   {
       fuel
       {
           type UO2;

           conductivity
           {
               type UO2Nfir;
           }

           heatCapacity
           {
               type UO2Cp;
           }
       }
   }

In this structure:

* the name of the subdictionary corresponds to the **material property**,
* the ``type`` entry selects the **runtime model** used to evaluate it.

This mechanism allows different models to be used for each property
independently.

-------------------------------------------------------------------------------

.. rubric:: Constant property shortcut

For materials with constant properties, OFFBEAT provides a simplified
syntax that avoids defining a full property subdictionary.

.. code-block:: cpp

   materials
   {
       cladding
       {
           type constant;

           k k             [1 -1 -3 -1 0 0 0] 16;
           Cp Cp           [0 2 -2 -1 0 0 0]  330;
           rho rho         [1 -3 0 0 0 0 0]   6500;
       }
   }

In this case:

* the property value is specified directly,
* no additional property model needs to be declared.

If a property requires a non-constant model, the dedicated
subdictionary form can still be used instead.

-------------------------------------------------------------------------------

.. rubric:: Separation of modeling layers

Within each material definition, OFFBEAT distinguishes between several
conceptual modeling layers.

1. Material properties
    Intrinsic properties such as density, thermal conductivity, or heat capacity.

2. Behavioral models
    Models describing material evolution driven by external conditions (e.g. swelling, densification, relocation).

3. Constitutive behavior
    Mechanical laws relating stress and strain (e.g. elasticity, plasticity, creep).

4. Damage models
    Models describing degradation and loss of load-carrying capacity (e.g. cracking, softening).

These layers are configured independently but act on the same material.

-------------------------------------------------------------------------------

.. rubric:: Writing material properties

For debugging or analysis purposes, the material properties computed by
the selected models can be written to disk.

This is enabled using the keyword ``writeMaterialProperties`` inside the
``materials`` dictionary.

.. code-block:: cpp

   materials
   {
       writeMaterialProperties on;

       fuel
       {
           type UO2;

           conductivity
           {
               type UO2Nfir;
           }
       }
   }

When this option is enabled, OFFBEAT writes the evaluated material
properties (e.g. conductivity, density, heat capacity, elastic
constants) as fields in the case directory. This allows users to inspect
the effective properties used by the solver and to verify correlations
or input settings.

The option is primarily intended for **verification, debugging, or
post-processing of material models** and is typically disabled in
production simulations.

-------------------------------------------------------------------------------

.. rubric:: Example

.. code-block:: cpp

   materials
   {
       fuel
       {
           type        UO2;
           enrichment  0.04;

           conductivity
           {
               type UO2Nfir;
           }

           swelling
           {
               type UO2Frapcon;
           }

           densification
           {
               type UO2Frapcon;
           }

           relocation
           {
               type UO2Frapcon;
           }

           rheology
           {
               type elasticity;
           }

           damage
           {
               type isotropicCracking;
           }
       }

       cladding
       {
           type Zircaloy;

           rho [1 -3 0 0 0 0 0] 6500;

           conductivity
           {
               type constant;
               value 16;
           }

           rheology
           {
               type elasticPlastic;
           }

           creep
           {
               type zircaloyCreep;
           }
       }
   }

In this example:

* ``fuel`` and ``cladding`` are ``cellZone`` names defined in the mesh,
* each zone is assigned exactly one material,
* properties, behavior, constitutive laws, and damage models are
  configured independently.

-------------------------------------------------------------------------------

The following pages describe the available material models.

.. toctree::
   :maxdepth: 1
   :caption: Contents

   materials/index
   properties/index
   behavior/index
   constitutive_behavior/index
   damage/index