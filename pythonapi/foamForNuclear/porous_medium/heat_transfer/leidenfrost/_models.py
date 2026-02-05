from foamForNuclear.common import OffbeatDict
from foamForNuclear._attrs_tools import offbeat_define
from typing import ClassVar
from attrs import field, validators as v

@offbeat_define
class LeidenfrostModel(OffbeatDict):
    TYPE: ClassVar[str] = "none"


@offbeat_define
class GroeneveldStewart(LeidenfrostModel):
    TYPE: ClassVar[str] = "GroeneveldStewart"
    criticalPressure: int | float


