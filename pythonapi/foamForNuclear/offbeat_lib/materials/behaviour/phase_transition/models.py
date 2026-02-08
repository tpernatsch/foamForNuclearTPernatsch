# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import ffn_define
from foamForNuclear.common import FoamForNuclearDict

@ffn_define
class PhaseTransitionModel(FoamForNuclearDict):
    """
    Base class for phase transition models.

    The parent class `type` is `none`. If selected, phase transition is not
    modeled, but the field `betaFraction` is created and registered.

    If a `betaFraction` file is present in the starting time folder, it is read.
    Otherwise, `betaFraction` is initialized to 0 with default boundary conditions.
    """
    TYPE: ClassVar[str] = 'none'

PhaseTransition = PhaseTransitionModel  # alias

@ffn_define
class ZircaloyDynamic(PhaseTransition):
    """
    Dynamic beta-phase transition model for Zircaloy-based claddings, derived from
    Massih (2009): "Transformation kinetics of zirconium alloys under non-isothermal
    conditions", Journal of Nuclear Materials 384 (2009) 330–335.

    This model evolves the `betaFraction` field based on the local temperature and
    time history.
    """
    TYPE: ClassVar[str] = 'ZircaloyDynamic'
