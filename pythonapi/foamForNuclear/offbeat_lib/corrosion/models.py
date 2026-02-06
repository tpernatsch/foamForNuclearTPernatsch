from __future__ import annotations

from foamForNuclear._attrs_tools import ffn_define
from typing import ClassVar
from foamForNuclear.common import FoamForNuclearDict, OpenFOAMDict
from attrs import field

@ffn_define
class Corrosion(FoamForNuclearDict):
    TYPE: ClassVar[str] = "constant"


# TODO: must be improved to allow  per-dict IDE suggestions
@ffn_define
class ByPatch(Corrosion):
    TYPE: ClassVar[str] = "byPatch"
    patches: OpenFOAMDict = field(factory=OpenFOAMDict)