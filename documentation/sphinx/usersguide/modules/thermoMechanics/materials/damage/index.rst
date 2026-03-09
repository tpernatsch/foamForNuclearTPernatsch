Damage models
=============

Damage models are used to represent the degradation of material integrity
due to cracking or damage accumulation. They capture the macroscopic
effects of damage on the mechanical response without explicitly resolving
individual cracks.

Depending on the specific model, damage effects may be represented by:

* the introduction of an additional damage-related strain contribution,
* a modification of the stress–strain relationship,
* a reduction of material stiffness or strength through softening laws.

Damage models are defined **per material** inside the ``materials``
dictionary of the ``solverDict``.

This configuration applies both when using the standalone **OFFBEAT**
application and when using the **thermo-mechanics module within a
multi-physics GeN-Foam simulation**.

Damage models interact with the constitutive behavior and rheology but
remain conceptually distinct from intrinsic material properties and from
behavioral models.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following damage models:

.. toctree::
   :maxdepth: 1
   :caption: Damage models
   :glob:

   ../../../../../../cppapi/generated/offbeatLib/materials/materialModel/damageModel/*