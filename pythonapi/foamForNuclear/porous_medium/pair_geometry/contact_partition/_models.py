from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import FoamForNuclearDict, OpenFOAMDict
from foamForNuclear._attrs_tools import ffn_define
from typing import ClassVar
from attrs import field

from collections.abc import Mapping

@ffn_define
class ContactPartitionModel(FoamForNuclearDict):
    TYPE: ClassVar[str] = "none"


@ffn_define
class ByRegime(ContactPartitionModel):
    TYPE: ClassVar[str] = "byRegime"
    regimeMap: str
    regimes: list[dict] = field(   # or list[dict] if you want, but object is safer here
        factory=list,
        metadata={"ffn_internal": True},
    )


    def _update_fields(self) -> None:
        super()._update_fields()

        for reg in self.regimes:
            if reg is None:
                continue

            # --- Case 1: flat dict API ---
            if isinstance(reg, Mapping):
                name = reg.get("name")
                if not name:
                    raise ValueError("Each regime dict must define a 'name'")

                # everything except 'name' goes into the OpenFOAM sub-dict
                items = {k: v for k, v in reg.items() if k != "name"}
                if not items:
                    raise ValueError(f"Regime '{name}' has no entries")

                reg_of = OpenFOAMDict(name=name, items=items)

            # --- Case 2: OpenFOAMDict-like object ---
            else:
                name = getattr(reg, "name", None)
                if not name:
                    raise ValueError("Each regime must have a 'name'")
                reg_of = reg

            self._of[name] = reg_of

    # TODO: Check if same logic can be replicated with some form of checked-dict type
    def add_regime(
            self,
            regimeName: str,
            regimeDict: OpenFOAMDict
        ) -> None:
        check_type("regimeName", regimeName, str)
        check_type("regimeDict", regimeDict, OpenFOAMDict)
        self.regimes[regimeName] = regimeDict

    def add_constant_regime(
            self,
            regimeName: str,
            value: float
        ) -> None:
        check_type("regimeName", regimeName, str)
        check_type("value", value, (float, int))
        self.add_regime(regimeName, OpenFOAMDict({
            "type": "constant",
            "value": value
        }))


@ffn_define
class Constant(ContactPartitionModel):
    TYPE: ClassVar[str] = "constant"
    value: float