Material models
=========

In OFFBEAT and in the thermo-mechanics module used within multi-physics
GeN-Foam simulations, a **material** defines the physical characteristics
of the solid assigned to a given mesh region.

Materials are associated with mesh regions through ``cellZone`` objects
and act as containers for:

* material properties (thermal and mechanical),
* behavioral models (e.g. swelling, densification),
* constitutive behavior (e.g. elasticity, plasticity),
* optional damage and failure models.

Materials are selected using the ``type`` keyword inside each material
subdictionary of the ``materials`` section in the ``solverDict``.

Each material comes with a **default set of thermo-mechanical properties
and associated models**, which reflect typical correlations and
assumptions for that material (e.g. UO₂, Zircaloy, steels).

These defaults can be selectively overridden by the user through
explicit model choices in the ``materials`` dictionary.

-------------------------------------------------------------------------------

.. rubric:: Material selection

.. code-block:: cpp

   materials
   {
       fuel
       {
           type UO2;
           // ...
           // material properties, behavior, constitutive models, etc.
       }
   }

The material type defines the **base material class**
(e.g. ``UO2``, ``Zircaloy``, ``Steel1515Ti``), independently of the
``cellZone`` name.

-------------------------------------------------------------------------------

.. rubric:: Scope of this section

This section documents the **available material types** in OFFBEAT and
in the thermo-mechanics module of GeN-Foam, together with their intended
usage.

Specific details on the various models for material properties,
behavioral phenomena, constitutive behavior, damage, and failure are
documented in the corresponding subsections of the user guide.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following material
types:

.. toctree::
   :maxdepth: 1
   :caption: Material classes
   :glob:

   ../../../../../../cppapi/generated/offbeatLib/materials/materialModel/*
   