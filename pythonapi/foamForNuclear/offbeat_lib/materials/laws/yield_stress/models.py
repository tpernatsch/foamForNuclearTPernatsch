from foamForNuclear.common import OffbeatDict
from foamForNuclear._attrs_tools import offbeat_define, _to_List_float
from typing import ClassVar
from foamForNuclear.common import Table, List
from attrs import field


@offbeat_define
class YieldStress(OffbeatDict):
    TYPE: ClassVar[str] = "none"


@offbeat_define
class Hardening(YieldStress):
    TYPE: ClassVar[str] = "hardening"
    plasticStrainVsYieldStress: Table | None = None


@offbeat_define
class Constant(YieldStress):
    TYPE: ClassVar[str] = "constant"
    sigmaY: float | int = 1.0


@offbeat_define
class FRAPTRAN(YieldStress):
    TYPE: ClassVar[str] = "Fraptran"
    phi: float | int
    CW: float | int 