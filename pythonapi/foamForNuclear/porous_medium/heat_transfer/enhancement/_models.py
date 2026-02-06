from foamForNuclear.common import FoamForNuclearDict
from foamForNuclear._attrs_tools import ffn_define
from typing import ClassVar
from attrs import field, validators as v

@ffn_define
class FlowEnhancementFactor(FoamForNuclearDict):
    TYPE: ClassVar[str] = "none"


@ffn_define
class CobraTf(FlowEnhancementFactor):
    TYPE: ClassVar[str] = "COBRA-TF"


@ffn_define
class RezkallahSims(FlowEnhancementFactor):
    TYPE: ClassVar[str] = "RezkallahSims"
    exp: int | float

