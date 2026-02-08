import numpy as np
from foamForNuclear.common import FoamForNuclearDict
from foamForNuclear._attrs_tools import ffn_define
from typing import ClassVar


@ffn_define
class TwoPhaseDragMultiplierModel(FoamForNuclearDict):
    TYPE: ClassVar[str] = "none"


@ffn_define
class LottesFlinn(TwoPhaseDragMultiplierModel):
    TYPE: ClassVar[str] = "LottesFlinn"
    multiplierFluid: str


@ffn_define
class LockhartMartinelli(TwoPhaseDragMultiplierModel):
    TYPE: ClassVar[str] = "LockhartMartinelli"
    multiplierFluid: str


@ffn_define
class Kaiser88(TwoPhaseDragMultiplierModel):
    TYPE: ClassVar[str] = "Kaiser88"
    multiplierFluid: str