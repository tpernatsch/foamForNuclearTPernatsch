# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import ffn_define
from foamForNuclear.common import FoamForNuclearDict

@ffn_define
class HeatCapacityModel(FoamForNuclearDict):
    """
    Mother class for heatCapacity models.
    """
    TYPE: ClassVar[str] = 'heatCapacity'

HeatCapacity = HeatCapacityModel  # alias

@ffn_define
class Steel1515TiBanerjee(HeatCapacity):
    """
    Correlation for 15-15 Ti heat capacity from Banerjee et al. (2007).


    Options
    -------
    par1 : scalar
        Model parameter.
        (default: 431; required: False)

    par2 : scalar
        Model parameter.
        (default: 0.177; required: False)

    par3 : scalar
        Model parameter.
        (default: 8.72e-05; required: False)
    """
    TYPE: ClassVar[str] = 'Steel1515TiBanerjee'
    par1: float | int = 431.0
    par2: float | int = 0.177
    par3: float | int = 8.72e-05

@ffn_define
class Constant(HeatCapacity):
    """
    Model for constant heat capacity. The value is read from dictionary.


    Options
    -------
    value : scalar
        Constant heat capacity value.
        (required: True)
    """
    TYPE: ClassVar[str] = 'constant'
    value: float | int

@ffn_define
class UPuO2Fink(HeatCapacity):
    """
    Heat Capacity model for (U,Pu)O2 derived from Fink.
    Source : https://info.ornl.gov/sites/publications/Files/Pub57523.pdf


    Options
    -------
    C1 : scalar
        Model parameter.
        (default: 322.49; required: False)

    C2 : scalar
        Model parameter.
        (default: 0.014679; required: False)

    C3 : scalar
        Model parameter.
        (default: 0; required: False)

    Ea : scalar
        Model parameter.
        (default: 18531.7; required: False)

    theta : scalar
        Model parameter.
        (default: 587.41; required: False)
    """
    TYPE: ClassVar[str] = 'UPuO2Fink'
    C1: float | int = 322.49
    C2: float | int = 0.014679
    C3: float | int = 0.0
    Ea: float | int = 18531.7
    theta: float | int = 587.41

@ffn_define
class ZircaloyIaea(HeatCapacity):
    """
    Correlation for Zircaloy heat capacity from Iaea.
    <a href="https://www-pub.iaea.org/MTCD/publications/PDF/te_1496_web.pdf">Click here</a>
    for the related Iaea publication.


    Options
    -------
    par1 : scalar
        Model parameter.
        (default: 255.66; required: False)

    par2 : scalar
        Model parameter.
        (default: 0.1024; required: False)

    par3 : scalar
        Model parameter.
        (default: 597.1; required: False)

    par4 : scalar
        Model parameter.
        (default: 0.4088; required: False)

    par5 : scalar
        Model parameter.
        (default: 0.0001565; required: False)

    par6 : scalar
        Model parameter.
        (default: 1058.4; required: False)

    par7 : scalar
        Model parameter.
        (default: 1213.8; required: False)

    par8 : scalar
        Model parameter.
        (default: 2.0; required: False)

    par9 : scalar
        Model parameter.
        (default: 719.61; required: False)
    """
    TYPE: ClassVar[str] = 'ZircaloyIaea'
    par1: float | int = 255.66
    par2: float | int = 0.1024
    par3: float | int = 597.1
    par4: float | int = 0.4088
    par5: float | int = 0.0001565
    par6: float | int = 1058.4
    par7: float | int = 1213.8
    par8: float | int = 2.0
    par9: float | int = 719.61

@ffn_define
class UO2Matpro(HeatCapacity):
    """
    Heat Capacity model for UO2 derived from Matprov11.


    Options
    -------
    ED : scalar
        Model parameter.
        (default: '1.577e5'; required: False)

    K1 : scalar
        Model parameter.
        (default: 296.7; required: False)

    K2 : scalar
        Model parameter.
        (default: 0.0243; required: False)

    K3 : scalar
        Model parameter.
        (default: '8.745e7'; required: False)

    OM : scalar
        Model parameter.
        (default: 2.0; required: False)

    par1 : scalar
        Model parameter.
        (default: 2.0; required: False)

    par2 : scalar
        Model parameter.
        (default: 2.0; required: False)

    par3 : scalar
        Model parameter.
        (default: 2.0; required: False)

    phi : scalar
        Model parameter.
        (default: 535.285; required: False)
    """
    TYPE: ClassVar[str] = 'UO2Matpro'
    ED: float | int = 1.577e5
    K1: float | int = 296.7
    K2: float | int = 0.0243
    K3: float | int = 8.745e7
    OM: float | int = 2.0
    par1: float | int = 2.0
    par2: float | int = 2.0
    par3: float | int = 2.0
    phi: float | int = 535.285

@ffn_define
class UPuO2Matpro(HeatCapacity):
    """
    Heat Capacity model for (U,Pu)O2 derived from Matprov11.


    Options
    -------
    ED : scalar
        Model parameter.
        (default: '1.967e5'; required: False)

    K1 : scalar
        Model parameter.
        (default: 347.4; required: False)

    K2 : scalar
        Model parameter.
        (default: '3.95e4'; required: False)

    K3 : scalar
        Model parameter.
        (default: '3.86e7'; required: False)

    OM : scalar
        Model parameter.
        (default: 2.0; required: False)

    oxygenMetalRatio : scalar
        Optional override. Unless strictly necessary it is suggested to define the
        oxygenMetalRatio in the base/composition dict.
        (required: False)

    ratioOverMetal : unknown
        Model parameter.
        (required: True)

    theta : scalar
        Model parameter.
        (default: 571; required: False)
    """
    TYPE: ClassVar[str] = 'UPuO2Matpro'
    ED: float | int = 1.967e5
    K1: float | int = 347.4
    K2: float | int = 3.95e4
    K3: float | int = 3.86e7
    OM: float | int = 2.0
    oxygenMetalRatio: float | int | None = None
    ratioOverMetal: Unknown
    theta: float | int = 571.0

@ffn_define
class ZircaloyMatpro(HeatCapacity):
    """
    Correlation for Zircaloy heat capacity from Matpro.
    Link : https://www.nrc.gov/docs/ML1429/ML14296A063.pdf - page 60.
    """
    TYPE: ClassVar[str] = 'ZircaloyMatpro'

@ffn_define
class Molybdenum(HeatCapacity):
    """
    Heat capacity correlation model `Molybdenum`.


    Options
    -------
    par1 : scalar
        Model parameter.
        (default: 9.74e-06; required: False)

    par2 : scalar
        Model parameter.
        (default: 0.0537; required: False)

    par3 : scalar
        Model parameter.
        (default: 235; required: False)
    """
    TYPE: ClassVar[str] = 'Molybdenum'
    par1: float | int = 9.74e-06
    par2: float | int = 0.0537
    par3: float | int = 235.0

@ffn_define
class SiCSnead(HeatCapacity):
    """
    Correlation for SiC heat capacity from Snead et al. (2007).


    Options
    -------
    par1 : scalar
        Model parameter.
        (default: 925.65; required: False)

    par2 : scalar
        Model parameter.
        (default: 0.3772; required: False)

    par3 : scalar
        Model parameter.
        (default: -7.9259e-05; required: False)

    par4 : scalar
        Model parameter.
        (default: '-3.19446e7'; required: False)
    """
    TYPE: ClassVar[str] = 'SiCSnead'
    par1: float | int = 925.65
    par2: float | int = 0.3772
    par3: float | int = -7.9259e-05
    par4: float | int = -3.19446e7
