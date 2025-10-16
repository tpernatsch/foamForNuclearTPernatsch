from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type


class PrghTotalPressure(Patch):
    """
    This boundary condition provides static pressure condition for `p_rgh`, calculated as:

        `p_rgh = p - ρg.(h - hRef)`

        `p = p0 - 0.5 ρ |U|²`

    where:
    - `p_rgh` = Pseudo hydrostatic pressure [Pa]
    - `p`     = Static pressure [Pa]
    - `p0`    = Total pressure [Pa]
    - `h`     = Height in the opposite direction to gravity
    - `hRef`  = Reference height in the opposite direction to gravity
    - `ρ`     = Density
    - `g`     = Acceleration due to gravity [m/s²]

    Parameters
    ----------
    value : float | int
        Value at the patch set uniformly
    p0 : float | int
        Total pressure
    U : str
        Velocity field name (default `None` = `U`)
    phi : str
        Flux field name (default `None` = `phi`)
    rho : str
        Density field name (default `thermo:rho`)
    """
    def __init__(
            self,
            value: float,
            p0: float,
            U: str=None,
            phi: str=None,
            rho: str="thermo:rho"
        ):
        super().__init__(type="prghTotalPressure", value=value)

        self.p0 = p0
        self.U = U
        self.phi = phi
        self.rho = rho

    @property
    def p0(self):
        return self._p0

    @p0.setter
    def p0(self, p0) -> None:
        check_type("p0", p0, (float, int))
        self._p0 = p0
        self.__setitem__('p0', f"uniform {p0}")

    @property
    def U(self):
        return self._U

    @U.setter
    def U(self, U) -> None:
        check_type("U", U, str, none_ok=True)
        self._U = U
        if (U is not None):
            self.__setitem__('U', U)

    @property
    def phi(self):
        return self._phi

    @phi.setter
    def phi(self, phi) -> None:
        check_type("phi", phi, str, none_ok=True)
        self._phi = phi
        if (phi is not None):
            self.__setitem__('phi', phi)

    @property
    def rho(self):
        return self._rho

    @rho.setter
    def rho(self, rho) -> None:
        check_type("rho", rho, str)
        self._rho = rho
        self.__setitem__('rho', rho)
