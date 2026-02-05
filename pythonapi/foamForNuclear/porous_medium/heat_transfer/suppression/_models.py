from foamForNuclear.common import OffbeatDict
from foamForNuclear._attrs_tools import offbeat_define
from typing import ClassVar
from attrs import field, validators as v

@offbeat_define
class SuppressionFactorModel(OffbeatDict):
    TYPE: ClassVar[str] = "none"


@offbeat_define
class CobraTf(SuppressionFactorModel):
    TYPE: ClassVar[str] = "COBRA-TF"


@offbeat_define
class Chen(SuppressionFactorModel):
    TYPE: ClassVar[str] = "Chen"

