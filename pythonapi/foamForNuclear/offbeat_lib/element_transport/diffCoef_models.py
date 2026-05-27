from __future__ import annotations
from typing import ClassVar
from foamForNuclear._attrs_tools import ffn_define
from foamForNuclear.common import FoamForNuclearDict


@ffn_define
class DiffCoefModel(FoamForNuclearDict):
    TYPE: ClassVar[str] = "diffCoefModel"
    failureTime: float | int | None = None


@ffn_define
class DoubleArrhenius(DiffCoefModel):
    TYPE: ClassVar[str] = "doubleArrhenius"
    d1: float | int
    q1: float | int
    d2: float | int
    q2: float | int


@ffn_define
class SingleArrhenius(DiffCoefModel):
    TYPE: ClassVar[str] = "singleArrhenius"
    d1: float | int
    q1: float | int


@ffn_define
class PiecewiseArrhenius(DiffCoefModel):
    TYPE: ClassVar[str] = "piecewiseArrhenius"
    d1: float | int
    q1: float | int
    d2: float | int
    q2: float | int
    Tswitch: float | int


@ffn_define
class PyCArrhenius(DiffCoefModel):
    TYPE: ClassVar[str] = "PyCArrhenius"


@ffn_define
class SiCArrhenius(DiffCoefModel):
    TYPE: ClassVar[str] = "SiCArrhenius"


@ffn_define
class UO2Arrhenius(DiffCoefModel):
    TYPE: ClassVar[str] = "UO2Arrhenius"


@ffn_define
class GraphiteArrhenius(DiffCoefModel):
    TYPE: ClassVar[str] = "graphiteArrhenius"
