from foamForNuclear.common import FoamForNuclearDict
from foamForNuclear._attrs_tools import ffn_define, _to_List_float
from typing import ClassVar
from foamForNuclear.common import Table, List
from attrs import field


@ffn_define
class YieldStress(FoamForNuclearDict):
    TYPE: ClassVar[str] = "none"


@ffn_define
class Hardening(YieldStress):
    TYPE: ClassVar[str] = "hardening"
    plasticStrainVsYieldStress: Table | None = None


@ffn_define
class Constant(YieldStress):
    TYPE: ClassVar[str] = "constant"
    sigmaY: float | int = 1.0


@ffn_define
class FRAPTRAN(YieldStress):
    TYPE: ClassVar[str] = "Fraptran"
    phi: float | int
    CW: float | int 