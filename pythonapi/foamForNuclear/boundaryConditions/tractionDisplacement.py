from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import Table, Vector
from foamForNuclear.timeProfile import OffbeatTimeProfile


class TractionDisplacement(Patch):
    """
    The `tractionDisplacement` fvPatchField can be selected in the patch subdictionary
    inside the `boundaryField` subdictionary of the displacement field.

    Parameters
    ----------
    tractionList : OffbeatTimeProfile
        A time-dependent traction table. Each entry consists of a time step and
        a corresponding traction vector value in Pa.
    traction : float
        Fixed traction vector in Pa. Used only if `tractionList` is absent.
    pressureList : OffbeatTimeProfile
        A time-dependent pressure table. Each entry consists of a time step and
        a corresponding pressure value in Pa.
    pressure : float
        Fixed pressure in Pa. Used only if `pressureList` is absent.
    planeStrain : bool
        Activates the plane strain approximation for the normal stress at the
        boundary. When enabled, the normal component of the strain (e.g.
        epsilonZZ) is assumed constant at the boundary. Useful for long
        axisymmetric rods. It cannot be used in combination with `flatSurface`
        keyword. **Default: `false`.**
    flatSurface : bool
        Activates the flat surface approximation for the normal stress at the
        boundary. When enabled, the normal component of the displacement (e.g.
        Dz) is assumed constant at the boundary, i.e. the patch is forced to
        remain flat. Useful for instance for the simulation of a half pellet
        fragment. It cannot be used in combination with `planeStrain` keyword.
        **Default: `false`.**
    fixedSpring : bool
        Activates a fixed spring-dashpot system for additional stability.
        **Default: `false`.**
    fixedSpringModulus : float
        Spring modulus in N/m. Required when `fixedSpring` is set to `true`.
    dashpotModulus : float
        Dashpot modulus in N/m. Required when `fixedSpring` is set to `true`.
    relax : float
        Relaxation factor for gradient updates. **Default: `1.0`.**
    value : Vector
        Initial displacement value (not stress).
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
        if (tractionList is not None):
            check_type("tractionList", tractionList, OffbeatTimeProfile)
            self.__setitem__('tractionList', tractionList)
        self._tractionList = tractionList

    @property
    def traction(self):
        return self._traction

    @traction.setter
    def traction(self, traction) -> None:
        check_type("traction", traction, Vector, none_ok=True)
        self._traction = traction
        if (traction is not None):
            self.__setitem__('traction', f"uniform {traction}")

    @property
    def pressureList(self):
        return self._pressureList

    @pressureList.setter
    def pressureList(self, pressureList) -> None:
        if (pressureList is not None):
            check_type("pressureList", pressureList, OffbeatTimeProfile)
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
