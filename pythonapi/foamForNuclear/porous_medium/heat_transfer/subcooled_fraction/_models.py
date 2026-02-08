from foamForNuclear.common import FoamForNuclearDict
from foamForNuclear._attrs_tools import ffn_define
from typing import ClassVar
from attrs import field, validators as v


@ffn_define
class SubCooledBoilingFractionModel(FoamForNuclearDict):
    TYPE: ClassVar[str] = "none"


@ffn_define
class SahaZuber(SubCooledBoilingFractionModel):
    TYPE: ClassVar[str] = "SahaZuber"


