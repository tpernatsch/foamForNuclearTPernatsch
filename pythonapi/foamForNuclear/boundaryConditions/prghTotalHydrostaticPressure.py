from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type


class PrghTotalHydrostaticPressure(Patch):
    """
    This boundary condition provides static pressure condition for `p_rgh`,
    calculated as

        `p_rgh = ph_rgh - 0.5 rho |U|2`

    where:

        - `p_rgh` = Pseudo hydrostatic pressure [Pa]
        - `ph_rgh` = Hydrostatic pressure: rho g(h-href) [Pa]
        - `h`     = Height in the opposite direction to gravity
        - `hRef`  = Reference height in the opposite direction to gravity
        - `rho`   = Density
        - `g`     = Acceleration due to gravity [m/s²]

    Parameters
    ----------
    value : float | int
        Value at the patch set uniformly
    ph_rgh : str
        ph_rgh field name
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
            ph_rgh: str=None,
            U: str=None,
            phi: str=None,
            rho: str="thermo:rho"
        ):
        super().__init__(type="prghTotalHydrostaticPressure", value=value)

        self.ph_rgh = ph_rgh
        self.U = U
        self.phi = phi
        self.rho = rho

    @property
    def ph_rgh(self):
        return self._ph_rgh

    @ph_rgh.setter
    def ph_rgh(self, ph_rgh) -> None:
        check_type("ph_rgh", ph_rgh, str)
        self._ph_rgh = ph_rgh
        if (ph_rgh is not None):
            self.__setitem__('ph_rgh', ph_rgh)

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
