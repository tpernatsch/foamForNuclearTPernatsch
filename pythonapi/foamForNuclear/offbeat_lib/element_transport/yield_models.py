from __future__ import annotations
from typing import ClassVar
from foamForNuclear._attrs_tools import ffn_define
from foamForNuclear.common import FoamForNuclearDict


@ffn_define
class YieldModel(FoamForNuclearDict):
    TYPE: ClassVar[str] = "none"


@ffn_define
class Constant(YieldModel):
    TYPE: ClassVar[str] = "constant"
    yield_: float | int

    def _update_fields(self) -> None:
        super()._update_fields()
        if "yield_" in self._of:
            self._of["yield"] = self._of.pop("yield_")


@ffn_define
class EnrichmentDependent(YieldModel):
    TYPE: ClassVar[str] = "enrichmentDependent"
    yieldLEU: float | int
    yieldHEU: float | int
    enrichmentThreshold: float | int = 0.175


@ffn_define
class BurnupDependent(YieldModel):
    TYPE: ClassVar[str] = "burnupDependent"
    yieldCoeff: float | int
    burnupExponent: float | int
    burnupConversion: float | int = 8.232


@ffn_define
class UO2(YieldModel):
    """Smart default: benchmark-standard yields for Cs/Ag/Kr/Sr in UO2 kernel."""
    TYPE: ClassVar[str] = "UO2"
