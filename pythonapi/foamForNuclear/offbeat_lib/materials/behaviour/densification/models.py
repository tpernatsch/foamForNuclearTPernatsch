# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import offbeat_define
from foamForNuclear.common import OffbeatDict

@offbeat_define
class DensificationModel(OffbeatDict):
    """
    Base class for densification models.

    The base model has type `none`. When selected, densification is not modeled,
    but the field `epsilonDensification` is still created and registered.

    If an `epsilonDensification` field file is provided in the starting time folder,
    it is read and used; otherwise the field is initialized to zero.
    """
    TYPE: ClassVar[str] = 'none'

Densification = DensificationModel  # alias

@offbeat_define
class Empirical(Densification):
    """
    Empirical densification model.

    The fuel density tends toward a final densification state as an exponential
    function of burnup. The user provides:
    - the final density increase due to densification
    - a time constant controlling how fast densification evolves with burnup


    Options
    -------
    burnupName : word
        Name of the burnup field to read from the mesh registry.
        (default: 'Bu'; required: False)

    densityFraction : scalar
        Initial density fraction (relative density) used by the model.
        (required: False)

    densificationDensityChange : scalar
        Final density increase due to densification in percent (model parameter).
        (required: True)

    densificationTimeConstant : scalar
        Time constant controlling densification evolution with burnup in MWd/t
        (model parameter).
        (required: True)
    """
    TYPE: ClassVar[str] = 'empirical'
    burnupName: str = 'Bu'
    densityFraction: float | int | None = None
    densificationDensityChange: float | int
    densificationTimeConstant: float | int

@offbeat_define
class UO2Frapcon(Densification):
    """
    Densification model derived from FRAPCON correlations for UO$_2$ fuel.


    Options
    -------
    burnupName : word
        Name of the burnup field to read from the mesh registry.
        (default: 'Bu'; required: False)

    densityFraction : scalar
        Initial density fraction (relative density) used by the model.
        (default: 0.95; required: False)

    Tsintering : scalar
        Sintering temperature used by the correlation [K].
        (default: 1800; required: False)

    resinteringDensityChange : scalar
        Density change associated with re-sintering (model parameter).
        (required: True)

    par1 : scalar
        Correlation parameter.
        (default: 22.2; required: False)

    par2 : scalar
        Correlation parameter.
        (default: 1453; required: False)

    par3 : scalar
        Correlation parameter.
        (default: 66.6; required: False)

    par4 : scalar
        Correlation parameter.
        (default: 2.0; required: False)

    par5 : scalar
        Correlation parameter.
        (default: 35.0; required: False)

    par6 : scalar
        Correlation parameter.
        (default: 3.0; required: False)
    """
    TYPE: ClassVar[str] = 'UO2Frapcon'
    burnupName: str = 'Bu'
    densityFraction: float | int = 0.95
    Tsintering: float | int = 1800.0
    resinteringDensityChange: float | int
    par1: float | int = 22.2
    par2: float | int = 1453.0
    par3: float | int = 66.6
    par4: float | int = 2.0
    par5: float | int = 35.0
    par6: float | int = 3.0
