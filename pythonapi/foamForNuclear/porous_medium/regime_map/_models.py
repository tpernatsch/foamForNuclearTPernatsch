from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import OffbeatDict, List, OpenFOAMDict
from foamForNuclear._attrs_tools import offbeat_define
from typing import ClassVar
from attrs import field, validators as v


_INTERPOLATION_MODE_TYPES = {"linear", "quadratic"}


@offbeat_define
class RegimeMapModel(OffbeatDict):
    """
    Base class for regime map model.
    """
    TYPE: ClassVar[str] = "none"


@offbeat_define
class OneParameter(RegimeMapModel):
    """
    Regime map model dependent on one parameter.
    """
    TYPE: ClassVar[str] = "oneParameter"
    parameter: str
    regimeBounds: dict[str, list] = field(factory=OpenFOAMDict)
    interpolationMode: str = field(
        default="linear",
        validator=v.and_(
            v.instance_of(str),
            v.in_(_INTERPOLATION_MODE_TYPES),
        ),
    )

    def add_regime(self, name: str, minValue: float, maxValue: float):
        check_type("name", name, str)
        check_type("minValue", minValue, (float, int))
        check_type("maxValue", maxValue, (float, int))

        self.regimeBounds[name] = List([minValue, maxValue])
