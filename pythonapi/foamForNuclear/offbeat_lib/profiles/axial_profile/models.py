from __future__ import annotations

from foamForNuclear._attrs_tools import ffn_define, _to_List_float
from typing import ClassVar
from foamForNuclear.common import FoamForNuclearDict, List, Table
from attrs import field

@ffn_define
class AxialProfile(FoamForNuclearDict):
    TYPE: ClassVar[str] = "none"


@ffn_define
class Flat(AxialProfile):
    TYPE: ClassVar[str] = "flat"


@ffn_define
class TimeDependentTabulated(AxialProfile):
    TYPE: ClassVar[str] = "timeDependentTabulated"
    timePoints: list | List = field(factory=list, converter=_to_List_float)
    axialLocations: list | List = field(factory=list, converter=_to_List_float)
    data: list | Table = field(converter=lambda v: Table(items=v, type=""))
    #to be improved the way we set default data, should be similar to timepoints
    axialInterpolationMethod: str = "linear"
    burnupInterpolationMethod: str = "linear"