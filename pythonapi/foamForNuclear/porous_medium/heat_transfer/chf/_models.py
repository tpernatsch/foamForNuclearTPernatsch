from foamForNuclear.common import OffbeatDict
from foamForNuclear._attrs_tools import offbeat_define
from typing import ClassVar
from attrs import field, validators as v


@offbeat_define
class CriticalHeatFluxModel(OffbeatDict):
    TYPE: ClassVar[str] = "none"


@offbeat_define
class ConstantCHF(CriticalHeatFluxModel):
    TYPE: ClassVar[str] = "constantCHF"
    value: int | float