from foamForNuclear.common import OffbeatDict, OpenFOAMDict
from foamForNuclear.common import OpenFOAMDict
from foamForNuclear._attrs_tools import offbeat_define
from typing import ClassVar


@offbeat_define
class InterfacialAreaDensityModel(OffbeatDict):
    TYPE: ClassVar[str] = "none"


@offbeat_define
class Spherical(InterfacialAreaDensityModel):
    TYPE: ClassVar[str] = "spherical"