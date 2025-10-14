from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_positive, check_type
from foamForNuclear.common import List, OpenFOAMDict, Table


class CoolantPressure(Patch):
    """
    Re-implementation of the tractionDisplacement mother class that allows the user
    to deal only with "coolantPressure" instead than the more confusion "traction" and
    "pressure".
    """
    def __init__(
            self,
            value: float,
            relax: float=1,
            coolantPressure: float=None,
            coolantPressureList: Table=None,
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
        check_type("coolantPressureList", coolantPressureList, Table, none_ok=True)
        self._coolantPressureList = coolantPressureList


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
