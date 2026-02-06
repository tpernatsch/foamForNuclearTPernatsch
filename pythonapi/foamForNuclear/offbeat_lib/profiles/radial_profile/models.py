from __future__ import annotations

from foamForNuclear._attrs_tools import ffn_define, _to_List_float
from typing import ClassVar
from foamForNuclear.common import FoamForNuclearDict, List, Table
from attrs import field

@ffn_define
class RadialProfile(FoamForNuclearDict):
    TYPE: ClassVar[str] = "none"


@ffn_define
class Flat(RadialProfile):
    TYPE: ClassVar[str] = "flat"


@ffn_define
class FromBurnup(RadialProfile):
    TYPE: ClassVar[str] = "fromBurnup"