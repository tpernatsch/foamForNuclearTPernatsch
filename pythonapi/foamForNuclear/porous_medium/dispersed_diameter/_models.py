from __future__ import annotations

from foamForNuclear.common import OffbeatDict
from foamForNuclear._attrs_tools import offbeat_define
from typing import ClassVar


@offbeat_define
class DispersedDiameterModel(OffbeatDict):
    TYPE: ClassVar[str] = "none"


@offbeat_define
class Constant(DispersedDiameterModel):
    TYPE: ClassVar[str] = "constant"
    value: int | float = 0.0