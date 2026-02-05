from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.boundaryConditions.tractionDisplacement import TractionDisplacement
from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import Table, Vector, OpenFOAMDict
from foamForNuclear.timeProfile import OffbeatTimeProfile


class TopCladRingPressure(TractionDisplacement):
    """
    Boundary condition that applies a **ring pressure** on the top of the cladding.

    The `topCladRingPressure` fvPatchField is selected in the patch subdictionary
    inside the `boundaryField` of the displacement field. It is intended to model
    the action of a pressure acting on an annular surface defined by an inner and
    outer radius (e.g. coolant or plenum pressure acting on the top cladding ring).

    The pressure can be prescribed either as a constant value or as a time-dependent
    table. For user convenience, the aliases `coolantPressure` and
    `coolantPressureList` can be used interchangeably with `pressure` and
    `pressureList`.

    Options
    -------
    value : list | Vector
        Initial displacement value on the patch.
        (required: True)

    topCapOuterRadius : scalar
        Outer radius of the annular surface where the pressure is applied.
        (required: True)

    topCapInnerRadius : scalar
        Inner radius of the annular surface where the pressure is applied.
        (required: True)

    tractionList : OffbeatTimeProfile
        Time-dependent traction table (Pa), inherited from
        `tractionDisplacement`.
        (required: False)

    traction : Vector
        Fixed traction vector (Pa). Used only if `tractionList` is absent.
        (required: False)

    pressureList : list[scalar] or Table
        Time-dependent pressure applied on the ring (Pa).
        Can also be provided via `coolantPressureList`.
        (required: False)

    pressure : scalar
        Fixed pressure applied on the ring (Pa).
        Can also be provided via `coolantPressure`.
        (required: False)

    coolantPressureList : list[scalar] or Table
        Alias for `pressureList`.
        (required: False)

    coolantPressure : scalar
        Alias for `pressure`.
        (required: False)

    outOfBounds : word
        Behaviour when querying a time-dependent pressure outside the provided
        time range (e.g. `clamp`).
        (default: clamp; required: False)

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
        topCapOuterRadius: float | int,
        topCapInnerRadius: float | int,
        tractionList: OffbeatTimeProfile = None,
        traction: Vector = None,
        pressureList: list[float] | Table=None,
        pressure: float = None,
        coolantPressureList: list[float] | Table=None,
        coolantPressure: float = None,
        outOfBounds: str='clamp',
        planeStrain: bool = False,
        flatSurface: bool = False,
        fixedSpring: bool = False,
        fixedSpringModulus: float = None,
        dashpotModulus: float = None,
        relax: float = 1,
    ):
        # ---- keyword alias resolution (coolantPressure* -> pressure*) ----
        if pressureList is not None and coolantPressureList is not None:
            raise TypeError("Use only one of pressureList or coolantPressureList")
        if pressure is not None and coolantPressure is not None:
            raise TypeError("Use only one of pressure or coolantPressure")

        if pressureList is None and coolantPressureList is not None:
            pressureList = coolantPressureList
        if pressure is None and coolantPressure is not None:
            pressure = coolantPressure

        super().__init__(
            value=value,
            tractionList=tractionList,
            traction=traction,
            pressureList=pressureList,
            pressure=pressure,
            planeStrain=planeStrain,
            flatSurface=flatSurface,
            fixedSpring=fixedSpring,
            fixedSpringModulus=fixedSpringModulus,
            dashpotModulus=dashpotModulus,
            relax=relax,
        )

        self.type = "topCladRingPressure"
        self.topCapOuterRadius = topCapOuterRadius
        self.topCapInnerRadius = topCapInnerRadius
        self.outOfBounds = outOfBounds


    def __repr__(self, depth: int=0):
        if self.pressureList is not None and self.coolantPressureList is None:
            self.coolantPressureList = self.pressureList
        if self.pressure is not None and self.coolantPressure is None:
            self.coolantPressure = self.pressure

        if (self.coolantPressure is not None and self.coolantPressureList is None):
            self.__setitem__("coolantPressure", f"uniform {self.coolantPressure}")

        elif (self.coolantPressure is None and self.coolantPressureList is not None):
            self.__setitem__("coolantPressureList", OpenFOAMDict({
                "coolantPressureList": self.coolantPressureList,
                "outOfBounds": self.outOfBounds
            }))

        else:
            msg = "In topCladRingPressure, provide either 'coolantPressure' or 'coolantPressureList'"
            raise ValueError(msg)

        self.__setitem__("relax", self.relax)

        return super().__repr__(depth)

    @property
    def topCapOuterRadius(self):
        return self._topCapOuterRadius

    @topCapOuterRadius.setter
    def topCapOuterRadius(self, topCapOuterRadius) -> None:
        check_type("topCapOuterRadius", topCapOuterRadius, (int, float))
        self._topCapOuterRadius = topCapOuterRadius
        self.__setitem__("topCapOuterRadius", topCapOuterRadius)

    @property
    def topCapInnerRadius(self):
        return self._topCapInnerRadius

    @topCapInnerRadius.setter
    def topCapInnerRadius(self, topCapInnerRadius) -> None:
        check_type("topCapInnerRadius", topCapInnerRadius, (int, float))
        self._topCapInnerRadius = topCapInnerRadius
        self.__setitem__("topCapInnerRadius", topCapInnerRadius)

    @property
    def outOfBounds(self):
        return self._outOfBounds

    @outOfBounds.setter
    def outOfBounds(self, outOfBounds) -> None:
        check_type("outOfBounds", outOfBounds, str)
        self._outOfBounds = outOfBounds

    # ---- keep canonical names (pressureList/pressure) ----
    @property
    def pressureList(self):
        return self._pressureList

    @pressureList.setter
    def pressureList(self, pressureList) -> None:
        if pressureList is not None:
            check_type("pressureList", pressureList, (list, Table))
            self._pressureList = Table(pressureList)
        self._pressureList = pressureList

    @property
    def pressure(self):
        return self._pressure

    @pressure.setter
    def pressure(self, pressure) -> None:
        check_type("pressure", pressure, (int, float), none_ok=True)
        self._pressure = pressure

    # ---- aliases (coolantPressure* <-> pressure*) ----
    @property
    def coolantPressureList(self):
        return Table(self.pressureList)

    @coolantPressureList.setter
    def coolantPressureList(self, v) -> None:
        self.pressureList = Table(v)

    @property
    def coolantPressure(self):
        return self.pressure

    @coolantPressure.setter
    def coolantPressure(self, v) -> None:
        self.pressure = v
