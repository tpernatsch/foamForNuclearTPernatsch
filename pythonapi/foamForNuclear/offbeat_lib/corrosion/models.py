from __future__ import annotations

from foamForNuclear._attrs_tools import offbeat_define
from typing import ClassVar
from foamForNuclear.common import OffbeatDict, OpenFOAMDict
from attrs import field

@offbeat_define
class Corrosion(OffbeatDict):
    TYPE: ClassVar[str] = "constant"


# TODO: must be improved to allow  per-dict IDE suggestions
@offbeat_define
class ByPatch(Corrosion):
    TYPE: ClassVar[str] = "byPatch"
    patches: OpenFOAMDict = field(factory=OpenFOAMDict)