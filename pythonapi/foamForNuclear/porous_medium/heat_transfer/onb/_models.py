from foamForNuclear.common import OffbeatDict
from foamForNuclear._attrs_tools import offbeat_define
from typing import ClassVar
from attrs import field, validators as v

@offbeat_define
class NucleateBoilingOnsetModel(OffbeatDict):
    TYPE: ClassVar[str] = "none"


@offbeat_define
class Basu(NucleateBoilingOnsetModel):
    TYPE: ClassVar[str] = "Basu"
    surfaceTension: int | float
    contactAngle: int | float

