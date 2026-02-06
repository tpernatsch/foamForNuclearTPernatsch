from foamForNuclear.common import FoamForNuclearDict
from foamForNuclear._attrs_tools import ffn_define
from typing import ClassVar
from attrs import field, validators as v

@ffn_define
class AnnularFlowModel(FoamForNuclearDict):
    TYPE: ClassVar[str] = "none"


@ffn_define
class CachardLiquid(AnnularFlowModel):
    TYPE: ClassVar[str] = "CachardLiquid"
    wallEmissivity: int | float
    liquidEmissivity: int | float


@ffn_define
class CachardVapour(AnnularFlowModel):
    TYPE: ClassVar[str] = "CachardVapour"

