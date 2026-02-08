# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import ffn_define
from foamForNuclear.common import FoamForNuclearDict

@ffn_define
class DensityModel(FoamForNuclearDict):
    """
    Mother class for density models.
    """
    TYPE: ClassVar[str] = 'density'

Density = DensityModel  # alias

@ffn_define
class MolybdenumConstant(Density):
    """
    Class to set Molybdenum density to a constant input value.


    Options
    -------
    densityValue : scalar
        Option read by the model.
        (default: 10280; required: False)
    """
    TYPE: ClassVar[str] = 'MolybdenumConstant'
    densityValue: float | int = 10280.0

@ffn_define
class UO2Constant(Density):
    """
    Class to set UO2 density to a constant input value.


    Options
    -------
    densityFraction : scalar
        Fraction of theoretical density.
        Optional override. If not provided here, the value is read from the
        material base dictionary (preferred).
        (required: False)

    theoreticalDensity : scalar
        Theoretical density.
        Optional override. If not provided here, the value is read from the
        material base dictionary (preferred).
        (required: False)
    """
    TYPE: ClassVar[str] = 'UO2Constant'
    densityFraction: float | int | None = None
    theoreticalDensity: float | int | None = None

@ffn_define
class UPuO2Constant(Density):
    """
    Class modelling the constant density of (U,Pu)O2 MOX fuel.


    Options
    -------
    densityFraction : scalar
        Fraction of theoretical density.
        Optional override. If not provided here, the value is read from the
        material base dictionary (preferred).
        (required: False)

    theoreticalDensity : scalar
        Theoretical density.
        Optional override. If not provided here, the value is read from the
        material base dictionary (preferred).
        (required: False)
    """
    TYPE: ClassVar[str] = 'UPuO2Constant'
    densityFraction: float | int | None = None
    theoreticalDensity: float | int | None = None

@ffn_define
class Constant(Density):
    """
    Model for constant density. The value is read from dictionary.


    Options
    -------
    value : scalar
        Density value in kg/m^3
        (required: True)
    """
    TYPE: ClassVar[str] = 'constant'
    value: float | int

@ffn_define
class ZircaloyIaea(Density):
    """
    Class modelling Zircaloy density from Iaea.


    Options
    -------
    par1 : scalar
        Option read by the model.
        (default: 6595.2; required: False)

    par2 : scalar
        Option read by the model.
        (default: 0.1477; required: False)

    par3 : scalar
        Option read by the model.
        (default: 6690.0; required: False)

    par4 : scalar
        Option read by the model.
        (default: 0.1855; required: False)
    """
    TYPE: ClassVar[str] = 'ZircaloyIaea'
    par1: float | int = 6595.2
    par2: float | int = 0.1477
    par3: float | int = 6690.0
    par4: float | int = 0.1855

@ffn_define
class Steel1515TiSchumann(Density):
    """
    Class modelling the density evolution of 15-15Ti cladding material through
    the Schumann (1970) correlation.


    Options
    -------
    par1 : scalar
        Option read by the model.
        (default: 7900.0; required: False)
    """
    TYPE: ClassVar[str] = 'Steel1515TiSchumann'
    par1: float | int = 7900.0
