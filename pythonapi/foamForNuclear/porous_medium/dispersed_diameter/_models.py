from __future__ import annotations

from foamForNuclear.common import FoamForNuclearDict
from foamForNuclear._attrs_tools import ffn_define
from typing import ClassVar


@ffn_define
class DispersedDiameterModel(FoamForNuclearDict):
    TYPE: ClassVar[str] = "none"


@ffn_define
class Constant(DispersedDiameterModel):
    TYPE: ClassVar[str] = "constant"
    value: int | float = 0.0