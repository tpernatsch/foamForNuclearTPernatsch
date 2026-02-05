from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_positive, check_type
from foamForNuclear.common import List, OpenFOAMDict, Table


class CoolantPressure(Patch):
    """
    Boundary condition that applies a coolant-side pressure to the solid, as a
    simplified re-implementation of OpenFOAM's `tractionDisplacement` interface.

    This wrapper lets the user specify only **coolantPressure**, avoiding the more
    confusing split between `traction` and `pressure` used by the original condition.

    The pressure can be provided either as a single constant value or as a list / table
    to define a time-dependent signal. When a time-dependent signal is used, values
    requested outside the provided range are handled according to `outOfBounds`.

    Options
    -------
    value : scalar
        Initial / placeholder boundary value for the displacement field (required by
        OpenFOAM patch-field syntax).
        (required: True)

    relax : scalar
        Relaxation factor applied when updating the imposed pressure.
        (default: 1; required: False)

    coolantPressure : scalar
        Constant coolant pressure to apply (Pa).
        Use this for a time-independent boundary condition.
        (required: False)

    coolantPressureList : list[scalar] or Table
        Time-dependent coolant pressure specification, either as a list of values or
        as a `Table` object.
        (required: False)

    outOfBounds : word
        Behaviour when querying the time-dependent coolant pressure outside the
        provided range (e.g. `clamp`).
        (default: clamp; required: False)
    """
    def __init__(
            self,
            value: float,
            relax: float=1,
            coolantPressure: float=None,
            coolantPressureList: list[float] | Table=None,
            outOfBounds: str='clamp'
        ):
        super().__init__(type="coolantPressure", value=value)

        self.relax = relax
        self.coolantPressure = coolantPressure
        self.coolantPressureList = coolantPressureList
        self.outOfBounds = outOfBounds


    def __repr__(self, depth: int=0):
        if (self.coolantPressure is not None and self.coolantPressureList is None):
            self.__setitem__("coolantPressure", f"uniform {self.coolantPressure}")

        elif (self.coolantPressure is None and self.coolantPressureList is not None):
            self.__setitem__("coolantPressureList", OpenFOAMDict({
                "coolantPressureList": self.coolantPressureList,
                "outOfBounds": self.outOfBounds
            }))

        else:
            msg = "In CoolantPressure, provide either 'coolantPressure' or 'coolantPressureList'"
            raise ValueError(msg)

        self.__setitem__("relax", self.relax)

        return super().__repr__(depth)


    @property
    def coolantPressure(self):
        return self._coolantPressure

    @coolantPressure.setter
    def coolantPressure(self, coolantPressure) -> None:
        check_type("coolantPressure", coolantPressure, (float, int), none_ok=True)
        self._coolantPressure = coolantPressure


    @property
    def coolantPressureList(self):
        return self._coolantPressureList

    @coolantPressureList.setter
    def coolantPressureList(self, coolantPressureList) -> None:
        check_type("coolantPressureList", coolantPressureList, (list, Table), none_ok=True)
        self._coolantPressureList = Table(coolantPressureList)


    @property
    def relax(self):
        return self._relax

    @relax.setter
    def relax(self, relax) -> None:
        check_type("relax", relax, (float, int))
        check_positive("relax", relax)
        self._relax = relax


    @property
    def outOfBounds(self):
        return self._outOfBounds

    @outOfBounds.setter
    def outOfBounds(self, outOfBounds) -> None:
        check_type("outOfBounds", outOfBounds, str)
        self._outOfBounds = outOfBounds
