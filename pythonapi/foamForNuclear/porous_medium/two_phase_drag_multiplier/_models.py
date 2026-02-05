import numpy as np
from foamForNuclear.common import OffbeatDict
from foamForNuclear._attrs_tools import offbeat_define
from typing import ClassVar


@offbeat_define
class TwoPhaseDragMultiplierModel(OffbeatDict):
    TYPE: ClassVar[str] = "none"


@offbeat_define
class LottesFlinn(TwoPhaseDragMultiplierModel):
    TYPE: ClassVar[str] = "LottesFlinn"
    multiplierFluid: str


@offbeat_define
class LockhartMartinelli(TwoPhaseDragMultiplierModel):
    TYPE: ClassVar[str] = "LockhartMartinelli"
    multiplierFluid: str


@offbeat_define
class Kaiser88(TwoPhaseDragMultiplierModel):
    TYPE: ClassVar[str] = "Kaiser88"
    multiplierFluid: str