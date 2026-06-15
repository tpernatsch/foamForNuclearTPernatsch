from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import Table, Vector
from foamForNuclear.timeProfile import OffbeatTimeProfile


def _to_time_profile(values):
    """Wrap a plain list in an OffbeatTimeProfile; return OffbeatTimeProfile unchanged.

    Each entry is (time, value) where value may be a list [x,y,z] (converted to
    Vector) or a scalar — as required by Table.
    """
    if isinstance(values, OffbeatTimeProfile):
        return values
    converted = [
        (t, Vector(v[0], v[1], v[2]) if isinstance(v, list) else v)
        for t, v in values
    ]
    return OffbeatTimeProfile('table', values=converted)


class TractionDisplacement(Patch):
    """
    Boundary condition that applies **traction and/or pressure** loading to the
    displacement field.

    The `tractionDisplacement` fvPatchField is selected in the patch subdictionary
    inside the `boundaryField` of the displacement field. It allows prescribing an
    external traction vector and/or a normal pressure on the patch, either as fixed
    values or as time-dependent profiles.

    Options
    -------
    value : list | Vector
        Initial displacement value on the patch.
        (required: True)

    tractionList : OffbeatTimeProfile
        Time-dependent traction table. Each entry provides a time and a traction
        vector value (Pa).
        (required: False)

    traction : Vector
        Fixed traction vector (Pa). Used only if `tractionList` is absent.
        (required: False)

    pressureList : OffbeatTimeProfile
        Time-dependent pressure table. Each entry provides a time and a pressure
        value (Pa).
        (required: False)

    pressure : scalar
        Fixed pressure (Pa). Used only if `pressureList` is absent.
        (required: False)

    planeStrain : bool
        Activate the plane strain approximation for the normal stress at the boundary.
        Cannot be used together with `flatSurface`.
        (default: False; required: False)

    flatSurface : bool
        Activate the flat surface approximation for the normal stress at the boundary.
        Cannot be used together with `planeStrain`.
        (default: False; required: False)

    fixedSpring : bool
        Enable a fixed spring-dashpot system for additional stability.
        (default: False; required: False)

    fixedSpringModulus : scalar
        Spring modulus in N/m (only if `fixedSpring` is true).
        (required: False)

    dashpotModulus : scalar
        Dashpot modulus in N/m (only if `fixedSpring` is true).
        (required: False)

    relax : scalar
        Relaxation factor for gradient updates.
        (default: 1; required: False)
    """

    def __init__(
            self,
            value,
            tractionList : OffbeatTimeProfile=None,
            traction: Vector=None,
            pressureList : OffbeatTimeProfile=None,
            pressure: float=None,
            planeStrain: bool=False,
            flatSurface: bool=False,
            fixedSpring: bool=False,
            fixedSpringModulus: float=None,
            dashpotModulus: float=None,
            relax: float=1
        ):
        super().__init__(type="tractionDisplacement", value=value)

        self.tractionList = tractionList
        self.traction = traction
        self.pressureList = pressureList
        self.pressure = pressure
        self.planeStrain = planeStrain
        self.flatSurface = flatSurface
        self.fixedSpring = fixedSpring
        self.fixedSpringModulus = fixedSpringModulus
        self.dashpotModulus = dashpotModulus
        self.relax = relax

    @property
    def tractionList(self):
        return self._tractionList

    @tractionList.setter
    def tractionList(self, tractionList) -> None:
        if tractionList is not None:
            check_type("tractionList", tractionList, (list, OffbeatTimeProfile))
            tractionList = _to_time_profile(tractionList)
            self.__setitem__('tractionList', tractionList)
        self._tractionList = tractionList

    @property
    def traction(self):
        return self._traction

    @traction.setter
    def traction(self, traction) -> None:
        check_type("traction", traction, (list, Vector), none_ok=True)
        if isinstance(traction, list):
            traction = Vector(traction[0], traction[1], traction[2])
        self._traction = traction
        if traction is not None:
            self.__setitem__('traction', f"uniform {traction}")

    @property
    def pressureList(self):
        return self._pressureList

    @pressureList.setter
    def pressureList(self, pressureList) -> None:
        if pressureList is not None:
            check_type("pressureList", pressureList, (list, OffbeatTimeProfile))
            pressureList = _to_time_profile(pressureList)
            self.__setitem__('pressureList', pressureList)
        self._pressureList = pressureList

    @property
    def pressure(self):
        return self._pressure

    @pressure.setter
    def pressure(self, pressure) -> None:
        check_type("pressure", pressure, (int, float), none_ok=True)
        self._pressure = pressure
        if (pressure is not None):
            self.__setitem__('pressure', f"uniform {pressure}")

    @property
    def planeStrain(self):
        return self._planeStrain

    @planeStrain.setter
    def planeStrain(self, planeStrain) -> None:
        check_type("planeStrain", planeStrain, bool)
        self._planeStrain = planeStrain
        self.__setitem__('planeStrain', planeStrain)

    @property
    def flatSurface(self):
        return self._flatSurface

    @flatSurface.setter
    def flatSurface(self, flatSurface) -> None:
        check_type("flatSurface", flatSurface, bool)
        self._flatSurface = flatSurface
        self.__setitem__('flatSurface', flatSurface)

    @property
    def fixedSpring(self):
        return self._fixedSpring

    @fixedSpring.setter
    def fixedSpring(self, fixedSpring) -> None:
        check_type("fixedSpring", fixedSpring, bool)
        self._fixedSpring = fixedSpring
        self.__setitem__('fixedSpring', fixedSpring)

    @property
    def relax(self):
        return self._relax

    @relax.setter
    def relax(self, relax) -> None:
        check_type("relax", relax, (int, float))
        self._relax = relax
        self.__setitem__('relax', relax)
