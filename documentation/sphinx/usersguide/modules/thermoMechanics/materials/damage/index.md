# <b>Damage Models</b>

Damage models are used to represent the degradation of material integrity due to
cracking or damage accumulation. They are intended to capture the macroscopic effects 
of damage on the mechanical response, without explicitly resolving individual cracks.

Depending on the specific model, damage effects may be represented by:
- the introduction of an additional damage-related strain contribution,
- a modification of the stress–strain relationship,
- a reduction of material stiffness or strength through softening laws.

Damage models are defined **per material** inside the `materials` dictionary of
the `solverDict` (located in the `constant` folder).

They interact with the constitutive behavior and rheology, but remain
conceptually distinct from intrinsic material properties and from behavioral
models.

## <b>Classes</b>

Currently, OFFBEAT supports the following damage models:

- [none](../../../classes/materials/materialModel/damageModel/damageModel.md)  
  no damage model applied.

- [isotropicCracking](../../../classes/materials/materialModel/damageModel/isotropicCracking.md)  
  isotropic damage model introducing stiffness degradation due to cracking.

- [mazars](../../../classes/materials/materialModel/damageModel/mazarsDamageModel.md)  
  Mazars-type damage model based on strain-driven damage evolution.
