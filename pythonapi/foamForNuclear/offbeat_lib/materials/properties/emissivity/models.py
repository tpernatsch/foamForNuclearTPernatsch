# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import offbeat_define
from foamForNuclear.common import OffbeatDict

@offbeat_define
class EmissivityModel(OffbeatDict):
    """
    Mother class for emissivity models.
    """
    TYPE: ClassVar[str] = 'emissivity'

Emissivity = EmissivityModel  # alias

@offbeat_define
class MolybdenumConstant(Emissivity):
    """
    Class to set Molybdenum emissivity to a constant input value.


    Options
    -------
    emissivityValue : scalar
        Constant emissivity value.
        (default: 0.2; required: False)
    """
    TYPE: ClassVar[str] = 'MolybdenumConstant'
    emissivityValue: float | int = 0.2

@offbeat_define
class ZircaloyConstant(Emissivity):
    """
    Class to set Zircaloy emissivity to a constant input value.


    Options
    -------
    emissivityValue : scalar
        Constant emissivity value.
        (default: 0.808642; required: False)
    """
    TYPE: ClassVar[str] = 'ZircaloyConstant'
    emissivityValue: float | int = 0.808642

@offbeat_define
class Constant(Emissivity):
    """
    Model for constant emissivity. The value is read from dictionary.


    Options
    -------
    value : scalar
        Constant emissivity value.
        (required: True)
    """
    TYPE: ClassVar[str] = 'constant'
    value: float | int

@offbeat_define
class UO2Relap(Emissivity):
    """
    UO2 emissivity model derived from Relap.


    Options
    -------
    par1 : scalar
        Correlation parameter.
        (default: 0.7856; required: False)

    par2 : scalar
        Correlation parameter.
        (default: 1.5263e-05; required: False)
    """
    TYPE: ClassVar[str] = 'UO2Relap'
    par1: float | int = 0.7856
    par2: float | int = 1.5263e-05
