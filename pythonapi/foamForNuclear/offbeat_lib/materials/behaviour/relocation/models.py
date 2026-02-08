# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import ffn_define
from foamForNuclear.common import FoamForNuclearDict

@ffn_define
class RelocationModel(FoamForNuclearDict):
    """
    Base class for fuel relocation models.

    The base model has type `none`. When selected, fuel relocation is not modeled,
    but the relocation-related fields are still created and registered.

    If a relocation field is provided in the initial time folder, it is read;
    otherwise it is initialized to zero.
    """
    TYPE: ClassVar[str] = 'none'

Relocation = RelocationModel  # alias

@ffn_define
class UO2Frapcon(Relocation):
    """
    Fuel relocation model based on the FRAPCON formulation.

    The model computes fuel relocation as a function of burnup and thermal state,
    and can include an optional recovery term.


    Options
    -------
    burnupName : word
        Name of the burnup field used by the relocation model.
        (default: 'Bu'; required: False)

    heatSourceName : word
        Name of the volumetric heat source field used by the model.
        (default: 'Q'; required: False)

    gapWidthName : word
        Name of the gap-width field used by the model.
        (default: 'gapWidth'; required: False)

    recoveryFraction : scalar
        Fraction controlling recovery of relocation (model parameter).
        (default: 0.0; required: False)

    relaxRecovery : scalar
        Relaxation factor applied to the recovery contribution.
        (default: 1.0; required: False)

    modifiedRelocationModel : bool
        Activates the modified FRAPCON relocation formulation.
        (required: False)

    outerPatch : str
        Name of outer fuel patch.
        (required: True)

    GapCold : scalar
        Cold gap width used by the model [m].
        (required: True)

    DiamCold : scalar
        Cold fuel diameter used by the model [m].
        (required: True)
    """
    TYPE: ClassVar[str] = 'UO2Frapcon'
    burnupName: str = 'Bu'
    heatSourceName: str = 'Q'
    gapWidthName: str = 'gapWidth'
    recoveryFraction: float | int = 0.0
    relaxRecovery: float | int = 1.0
    modifiedRelocationModel: bool | None = None
    outerPatch: str
    GapCold: float | int
    DiamCold: float | int
