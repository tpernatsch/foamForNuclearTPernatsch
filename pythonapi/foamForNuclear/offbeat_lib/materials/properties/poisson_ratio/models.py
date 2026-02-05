# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import offbeat_define
from foamForNuclear.common import OffbeatDict

@offbeat_define
class PoissonRatioModel(OffbeatDict):
    """
    Mother class for PoissonRatio models.
    """
    TYPE: ClassVar[str] = 'PoissonRatio'

PoissonRatio = PoissonRatioModel  # alias

@offbeat_define
class MolybdenumConstant(PoissonRatio):
    """
    Class to set Mo Poisson's Ratio to a constant input value.


    Options
    -------
    PoissonRatioValue : scalar
        Option read by the model.
        (default: 0.316; required: False)
    """
    TYPE: ClassVar[str] = 'MolybdenumConstant'
    PoissonRatioValue: float | int = 0.316

@offbeat_define
class UO2Constant(PoissonRatio):
    """
    Class to set UO2 Poisson's Ratio to a constant input value.


    Options
    -------
    PoissonRatioValue : scalar
        Option read by the model.
        (default: 0.316; required: False)
    """
    TYPE: ClassVar[str] = 'UO2Constant'
    PoissonRatioValue: float | int = 0.316

@offbeat_define
class UPuO2Constant(PoissonRatio):
    """
    Class modelling the constant Poisson's ratio of UPuO2 MOX fuel from
    Matprov11.


    Options
    -------
    PoissonRatioValue : scalar
        Option read by the model.
        (default: 0.276; required: False)
    """
    TYPE: ClassVar[str] = 'UPuO2Constant'
    PoissonRatioValue: float | int = 0.276

@offbeat_define
class ZircaloyConstant(PoissonRatio):
    """
    Class to set Zircaloy Poisson's Ratio to a constant input value.


    Options
    -------
    PoissonRatioValue : scalar
        Option read by the model.
        (default: 0.3; required: False)
    """
    TYPE: ClassVar[str] = 'ZircaloyConstant'
    PoissonRatioValue: float | int = 0.3

@offbeat_define
class Constant(PoissonRatio):
    """
    Model for constant Poisson Ratio. The value is read from dictionary.


    Options
    -------
    value : scalar
        Poisson ratio value in [-]
        (required: True)
    """
    TYPE: ClassVar[str] = 'constant'
    value: float | int

@offbeat_define
class ZircaloyMatpro(PoissonRatio):
    """
    Model for Zircaloy Poisson's Ratio from Matprov11.


    Options
    -------
    YoungModulusName : word
        Option read by the model.
        (default: 'E'; required: False)

    fastFluenceName : word
        Option read by the model.
        (default: 'fastFluence'; required: False)

    par1 : scalar
        Option read by the model.
        (default: '7.07e11'; required: False)

    par10 : scalar
        Option read by the model.
        (default: '1.66e7'; required: False)

    par2 : scalar
        Option read by the model.
        (default: '2.315e8'; required: False)

    par3 : scalar
        Option read by the model.
        (default: '2.6e10'; required: False)

    par4 : scalar
        Option read by the model.
        (default: 0.88; required: False)

    par5 : scalar
        Option read by the model.
        (default: 0.12; required: False)

    par6 : scalar
        Option read by the model.
        (default: '1e25'; required: False)

    par7 : scalar
        Option read by the model.
        (default: '4.04e10'; required: False)

    par8 : scalar
        Option read by the model.
        (default: '2.168e7'; required: False)

    par9 : scalar
        Option read by the model.
        (default: '3.49e10'; required: False)
    """
    TYPE: ClassVar[str] = 'ZircaloyMatpro'
    YoungModulusName: str = 'E'
    fastFluenceName: str = 'fastFluence'
    par1: float | int = 7.07e11
    par10: float | int = 1.66e7
    par2: float | int = 2.315e8
    par3: float | int = 2.6e10
    par4: float | int = 0.88
    par5: float | int = 0.12
    par6: float | int = 1e25
    par7: float | int = 4.04e10
    par8: float | int = 2.168e7
    par9: float | int = 3.49e10

@offbeat_define
class Steel1515TiTobbe(PoissonRatio):
    """
    Class modelling Poisson's ratio of 15-15 Ti cladding material based on Tobbe
    correlation (1975).


    Options
    -------
    par1 : scalar
        Option read by the model.
        (default: 0.277; required: False)

    par2 : scalar
        Option read by the model.
        (default: '6e-5'; required: False)
    """
    TYPE: ClassVar[str] = 'Steel1515TiTobbe'
    par1: float | int = 0.277
    par2: float | int = 6e-5
