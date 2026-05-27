from __future__ import annotations
from typing import ClassVar
from attrs import field
from foamForNuclear._attrs_tools import ffn_define
from foamForNuclear.common import FoamForNuclearDict, OpenFOAMDict, List
from foamForNuclear.offbeat_lib.misc import ElementTransportSolver


def _zone_dict(mapping: dict) -> OpenFOAMDict:
    d = OpenFOAMDict()
    for zone, val in mapping.items():
        if isinstance(val, dict):
            inner = OpenFOAMDict()
            for sp, m in val.items():
                inner[sp] = m
            d[zone] = inner
        else:
            d[zone] = val
    return d


@ffn_define
class FpDiffusion(ElementTransportSolver):
    TYPE: ClassVar[str] = "byList"
    fissionProducts: list = field(factory=list, metadata={"ffn_internal": True})
    diffusion: dict       = field(factory=dict, metadata={"ffn_internal": True})
    yield_: dict          = field(factory=dict, metadata={"ffn_internal": True})

    def _update_fields(self) -> None:
        self._of.clear()
        self._of["type"] = "byList"
        self._of["solvers"] = List(["FpDiffusion"])

        fp_opts = OpenFOAMDict()
        fp_opts["fissionProducts"] = List(self.fissionProducts)
        fp_opts["diffusion"] = _zone_dict(self.diffusion)
        if self.yield_:
            fp_opts["yield"] = _zone_dict(self.yield_)

        self._of["FpDiffusionOptions"] = fp_opts
