from foamForNuclear.common import FoamForNuclearDict, OpenFOAMDict
from foamForNuclear.common import OpenFOAMDict
from foamForNuclear._attrs_tools import ffn_define
from typing import ClassVar


@ffn_define
class InterfacialAreaDensityModel(FoamForNuclearDict):
    TYPE: ClassVar[str] = "none"


@ffn_define
class Spherical(InterfacialAreaDensityModel):
    TYPE: ClassVar[str] = "spherical"