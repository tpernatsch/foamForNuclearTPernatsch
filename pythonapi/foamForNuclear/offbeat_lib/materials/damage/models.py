# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import offbeat_define
from foamForNuclear.common import OffbeatDict

@offbeat_define
class DamageModel(OffbeatDict):
    """
    Base class for damage models. The base model has type `none`. When selected, damage is
    not modeled.
    """
    TYPE: ClassVar[str] = 'none'

Damage = DamageModel  # alias

@offbeat_define
class IsotropicCracking(Damage):
    """
    Isotropic cracking damage model.

    The model supports alternative isotropic cracking formulations selected via
    `isotropicCrackingType`.


    Options
    -------
    isotropicCrackingType : word
        Choice of isotropic cracking formulation. Either `Barani` or `JankusWeeks`
        (default: 'Barani'; required: False)

    nCracksMax : scalar
        Maximum number of cracks (from 0 to 12) for the Barani model.
        If the JankusWeeks model is selected, this keyword is interepreted as the
        fixed number of cracks.
        (default: 12; required: False)
    """
    TYPE: ClassVar[str] = 'isotropicCracking'
    isotropicCrackingType: str = 'Barani'
    nCracksMax: float | int = 12.0

@offbeat_define
class Mazars(Damage):
    """
    Mazars scalar damage model.


    Options
    -------
    excludeAxialCmpt : bool
        If true, excludes the axial stress/strain component from the equivalent
        strain used to drive damage.
        (default: False; required: False)

    relaxDamage : scalar
        Relaxation factor applied to the damage update.
        (default: 1.0; required: False)

    d_max : scalar
        Maximum allowed value for damage variable.
        (default: 0.99; required: False)

    temperatureFragmentationThreshold : scalar
        Temperature threshold above wwhich fragmentation starts.
        (default: 1500; required: False)

    equivalentStrainRateThreshold : bool
        Activate equivalent strain rate threshold criterion
        (default: False; required: False)

    strainRateThreshold : scalar
        Equivalent strain rate threshold (valid only if previous
        `equivalentStrainRateThreshold` is set to true)
        (default: 0.1; required: False)
    """
    TYPE: ClassVar[str] = 'mazars'
    excludeAxialCmpt: bool = False
    relaxDamage: float | int = 1.0
    d_max: float | int = 0.99
    temperatureFragmentationThreshold: float | int = 1500.0
    equivalentStrainRateThreshold: bool = False
    strainRateThreshold: float | int = 0.1
