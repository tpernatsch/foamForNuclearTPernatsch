# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import offbeat_define
from foamForNuclear.common import OffbeatDict

@offbeat_define
class YoungModulusModel(OffbeatDict):
    """
    Mother class for YoungModulus models.
    """
    TYPE: ClassVar[str] = 'YoungModulus'

YoungModulus = YoungModulusModel  # alias

@offbeat_define
class Constant(YoungModulus):
    """
    Model for constant Young modulus. The value is read from dictionary.


    Options
    -------
    value : scalar
        Constant Young Modulus in Pa.
        (required: True)
    """
    TYPE: ClassVar[str] = 'constant'
    value: float | int

@offbeat_define
class SteelD9Hofman(YoungModulus):
    """
    Class modelling the Young's modulus of D9 Steel cladding using the Hofman
    correlation.


    Options
    -------
    par1 : scalar
        Option read by the model.
        (default: '2.01e5'; required: False)

    par2 : scalar
        Option read by the model.
        (default: 79.29; required: False)
    """
    TYPE: ClassVar[str] = 'SteelD9Hofman'
    par1: float | int = 2.01e5
    par2: float | int = 79.29

@offbeat_define
class UO2Matpro(YoungModulus):
    """
    Class modelling Young Modulus of UO2 fuel from Matprov11.


    Options
    -------
    densityFraction : scalar
        Fraction of theoretical density.
        Optional override. If not provided here, the value is read from the
        material base dictionary (preferred).
        (required: False)

    par1 : scalar
        Option read by the model.
        (default: '2.334e11'; required: False)

    par2 : scalar
        Option read by the model.
        (default: 2.752; required: False)

    par3 : scalar
        Option read by the model.
        (default: 0.00010915; required: False)
    """
    TYPE: ClassVar[str] = 'UO2Matpro'
    densityFraction: float | int | None = None
    par1: float | int = 2.334e11
    par2: float | int = 2.752
    par3: float | int = 0.00010915

@offbeat_define
class UPuO2Matpro(YoungModulus):
    """
    Class modelling the Young's modulus of (U,Pu)O2 MOX fuel.


    Options
    -------
    densityFraction : scalar
        Fraction of theoretical density.
        Optional override. If not provided here, the value is read from the
        material base dictionary (preferred).
        (required: False)

    oxygenMetalRatio : scalar
        Oxygen to metal ratio.
        Optional override. If not provided here, the value is read from the
        material base dictionary (preferred).
        (required: False)

    par1 : scalar
        Option read by the model.
        (default: '2.334e11'; required: False)

    par2 : scalar
        Option read by the model.
        (default: 2.752; required: False)

    par3 : scalar
        Option read by the model.
        (default: 0.00010915; required: False)
    """
    TYPE: ClassVar[str] = 'UPuO2Matpro'
    densityFraction: float | int | None = None
    oxygenMetalRatio: float | int | None = None
    par1: float | int = 2.334e11
    par2: float | int = 2.752
    par3: float | int = 0.00010915

@offbeat_define
class ZircaloyMatpro(YoungModulus):
    """
    Class modelling Young Modulus of Zircaloy from Matprov11.


    Options
    -------
    fastFluenceName : word
        Option read by the model.
        (default: 'fastFluence'; required: False)

    par1 : scalar
        Option read by the model.
        (default: '6.61e11'; required: False)

    par10 : scalar
        Option read by the model.
        (default: '4.05e7'; required: False)

    par2 : scalar
        Option read by the model.
        (default: '5.912e8'; required: False)

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
        (default: '1.088e11'; required: False)

    par8 : scalar
        Option read by the model.
        (default: '5.475e7'; required: False)

    par9 : scalar
        Option read by the model.
        (default: '9.21e10'; required: False)
    """
    TYPE: ClassVar[str] = 'ZircaloyMatpro'
    fastFluenceName: str = 'fastFluence'
    par1: float | int = 6.61e11
    par10: float | int = 4.05e7
    par2: float | int = 5.912e8
    par3: float | int = 2.6e10
    par4: float | int = 0.88
    par5: float | int = 0.12
    par6: float | int = 1e25
    par7: float | int = 1.088e11
    par8: float | int = 5.475e7
    par9: float | int = 9.21e10

@offbeat_define
class Molybdenum(YoungModulus):
    """
    Class modelling Young Modulus of Molybdenum from Bison manual.


    Options
    -------
    par1 : scalar
        Option read by the model.
        (default: 334900000000.0; required: False)

    par2 : scalar
        Option read by the model.
        (default: 51010000.0; required: False)
    """
    TYPE: ClassVar[str] = 'Molybdenum'
    par1: float | int = 334900000000.0
    par2: float | int = 51010000.0

@offbeat_define
class BufferParfume(YoungModulus):
    """
    Class modelling Young Modulus of Bufffer from Parfume.


    Options
    -------
    densityName : word
        Option read by the model.
        (default: 'rho'; required: False)

    fastFluenceName : word
        Option read by the model.
        (default: 'fastFluence'; required: False)

    par1 : scalar
        Option read by the model.
        (default: 25.5; required: False)

    par2 : scalar
        Option read by the model.
        (default: 0.384; required: False)

    par3 : scalar
        Option read by the model.
        (default: 0.000324; required: False)

    par4 : scalar
        Option read by the model.
        (default: 0.23; required: False)

    par5 : scalar
        Option read by the model.
        (default: 0.00015; required: False)

    par6 : scalar
        Option read by the model.
        (default: 20.0; required: False)
    """
    TYPE: ClassVar[str] = 'BufferParfume'
    densityName: str = 'rho'
    fastFluenceName: str = 'fastFluence'
    par1: float | int = 25.5
    par2: float | int = 0.384
    par3: float | int = 0.000324
    par4: float | int = 0.23
    par5: float | int = 0.00015
    par6: float | int = 20.0

@offbeat_define
class PyCParfume(YoungModulus):
    """
    Class modelling Young Modulus of PyC from Parfume.


    Options
    -------
    densityName : word
        Option read by the model.
        (default: 'rho'; required: False)

    fastFluenceName : word
        Option read by the model.
        (default: 'fastFluence'; required: False)

    par1 : scalar
        Option read by the model.
        (default: 25.5; required: False)

    par10 : scalar
        Option read by the model.
        (default: 20.0; required: False)

    par11 : scalar
        Option read by the model.
        (default: 0.481; required: False)

    par12 : scalar
        Option read by the model.
        (default: 0.519; required: False)

    par2 : scalar
        Option read by the model.
        (default: 0.384; required: False)

    par3 : scalar
        Option read by the model.
        (default: 0.000324; required: False)

    par4 : scalar
        Option read by the model.
        (default: 1.463; required: False)

    par5 : scalar
        Option read by the model.
        (default: -0.463; required: False)

    par6 : scalar
        Option read by the model.
        (default: 2.985; required: False)

    par7 : scalar
        Option read by the model.
        (default: -0.0662; required: False)

    par8 : scalar
        Option read by the model.
        (default: 0.23; required: False)

    par9 : scalar
        Option read by the model.
        (default: 0.00015; required: False)
    """
    TYPE: ClassVar[str] = 'PyCParfume'
    densityName: str = 'rho'
    fastFluenceName: str = 'fastFluence'
    par1: float | int = 25.5
    par10: float | int = 20.0
    par11: float | int = 0.481
    par12: float | int = 0.519
    par2: float | int = 0.384
    par3: float | int = 0.000324
    par4: float | int = 1.463
    par5: float | int = -0.463
    par6: float | int = 2.985
    par7: float | int = -0.0662
    par8: float | int = 0.23
    par9: float | int = 0.00015

@offbeat_define
class SiCParfume(YoungModulus):
    """
    Class modelling Young Modulus of SiC with the Interpolation from Parfume.


    Options
    -------
    E1 : scalar
        Option read by the model.
        (default: 428.0; required: False)

    E2 : scalar
        Option read by the model.
        (default: 375.0; required: False)

    E3 : scalar
        Option read by the model.
        (default: 340.0; required: False)

    E4 : scalar
        Option read by the model.
        (default: 198.0; required: False)

    temp1 : scalar
        Option read by the model.
        (default: 25.0; required: False)

    temp2 : scalar
        Option read by the model.
        (default: 940.0; required: False)

    temp3 : scalar
        Option read by the model.
        (default: 1215.0; required: False)

    temp4 : scalar
        Option read by the model.
        (default: 1600.0; required: False)
    """
    TYPE: ClassVar[str] = 'SiCParfume'
    E1: float | int = 428.0
    E2: float | int = 375.0
    E3: float | int = 340.0
    E4: float | int = 198.0
    temp1: float | int = 25.0
    temp2: float | int = 940.0
    temp3: float | int = 1215.0
    temp4: float | int = 1600.0

@offbeat_define
class UPuO2SckCen(YoungModulus):
    """
    Correlation for Young's modulus of (U,Pu)O2 MOX fuel from SCK-CEN
    developed for TRANSURANUS.

    Source:
    https://re.public.polimi.it/retrieve/handle/11311/1172415/619939/Technical%20Report_INSPYRE_WP7-D7.2_%282020%29.pdf


    Options
    -------
    E1 : scalar
        Option read by the model.
        (default: 0.9796; required: False)

    E2 : scalar
        Option read by the model.
        (default: 3.1597e-05; required: False)

    E3 : scalar
        Option read by the model.
        (default: -7.6164e-08; required: False)

    EPuO2 : scalar
        Option read by the model.
        (default: 249.45; required: False)

    EUO2 : scalar
        Option read by the model.
        (default: 218.74; required: False)

    oxygenMetalRatio : scalar
        Oxygen to metal ratio.
        Optional override. If not provided here, the value is read from the
        material base dictionary (preferred).
        (required: False)
    """
    TYPE: ClassVar[str] = 'UPuO2SckCen'
    E1: float | int = 0.9796
    E2: float | int = 3.1597e-05
    E3: float | int = -7.6164e-08
    EPuO2: float | int = 249.45
    EUO2: float | int = 218.74
    oxygenMetalRatio: float | int | None = None

@offbeat_define
class SiCSnead(YoungModulus):
    """
    Class modelling Young Modulus of SiC from Snead et al. work:
    L. L. Snead, T. Nozawa, Y. Katoh, T.-S. Byun, S. Kondo, and D. A. Petti.
    Handbook of sic properties for fuel performance modeling.
    Journal of Nuclear Materials, 371:329–377, 2007.


    Options
    -------
    B : scalar
        Option read by the model.
        (default: '0.04e9'; required: False)

    C : scalar
        Option read by the model.
        (default: 0.0; required: False)

    E0 : scalar
        Option read by the model.
        (default: '460e9'; required: False)

    T0 : scalar
        Option read by the model.
        (default: 962.0; required: False)
    """
    TYPE: ClassVar[str] = 'SiCSnead'
    B: float | int = 0.04e9
    C: float | int = 0.0
    E0: float | int = 460e9
    T0: float | int = 962.0

@offbeat_define
class Steel1515TiTobbe(YoungModulus):
    """
    Class modelling the Young's modulus of 15-15 Ti using the Tobbe correlation
    (1975).


    Options
    -------
    par1 : scalar
        Option read by the model.
        (default: 202.7; required: False)

    par2 : scalar
        Option read by the model.
        (default: 0.08167; required: False)
    """
    TYPE: ClassVar[str] = 'Steel1515TiTobbe'
    par1: float | int = 202.7
    par2: float | int = 0.08167

@offbeat_define
class HastelloyNWatrous(YoungModulus):
    """
    Class modelling Young Modulus of HN from work of Watrous.


    Options
    -------
    par1 : scalar
        Option read by the model.
        (default: 9.944e-08; required: False)

    par2 : scalar
        Option read by the model.
        (default: 0.0001178; required: False)

    par3 : scalar
        Option read by the model.
        (default: 0.1033; required: False)

    par4 : scalar
        Option read by the model.
        (default: 220.9; required: False)
    """
    TYPE: ClassVar[str] = 'HastelloyNWatrous'
    par1: float | int = 9.944e-08
    par2: float | int = 0.0001178
    par3: float | int = 0.1033
    par4: float | int = 220.9
