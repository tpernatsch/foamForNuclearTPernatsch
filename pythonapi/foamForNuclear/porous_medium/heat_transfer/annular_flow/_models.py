from foamForNuclear.common import OffbeatDict
from foamForNuclear._attrs_tools import offbeat_define
from typing import ClassVar
from attrs import field, validators as v

@offbeat_define
class AnnularFlowModel(OffbeatDict):
    TYPE: ClassVar[str] = "none"


@offbeat_define
class CachardLiquid(AnnularFlowModel):
    TYPE: ClassVar[str] = "CachardLiquid"
    wallEmissivity: int | float
    liquidEmissivity: int | float


@offbeat_define
class CachardVapour(AnnularFlowModel):
    TYPE: ClassVar[str] = "CachardVapour"

