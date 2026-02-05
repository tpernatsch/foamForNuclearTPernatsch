from foamForNuclear.common import OffbeatDict
from foamForNuclear._attrs_tools import offbeat_define
from typing import ClassVar
from attrs import field, validators as v

@offbeat_define
class FlowEnhancementFactor(OffbeatDict):
    TYPE: ClassVar[str] = "none"


@offbeat_define
class CobraTf(FlowEnhancementFactor):
    TYPE: ClassVar[str] = "COBRA-TF"


@offbeat_define
class RezkallahSims(FlowEnhancementFactor):
    TYPE: ClassVar[str] = "RezkallahSims"
    exp: int | float

