from __future__ import annotations

from foamForNuclear._attrs_tools import offbeat_define, _to_List_float
from typing import ClassVar
from foamForNuclear.common import OffbeatDict, List, Table
from attrs import field

@offbeat_define
class RadialProfile(OffbeatDict):
    TYPE: ClassVar[str] = "none"


@offbeat_define
class Flat(RadialProfile):
    TYPE: ClassVar[str] = "flat"


@offbeat_define
class FromBurnup(RadialProfile):
    TYPE: ClassVar[str] = "fromBurnup"