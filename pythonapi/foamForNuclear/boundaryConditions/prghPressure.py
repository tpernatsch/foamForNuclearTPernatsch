from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type


class PrghPressure(Patch):
    """
    This boundary condition provides static pressure condition for p_rgh, calculated as:

        p_rgh = p - ρg.(h - hRef)

    where:
    - prgh  = Pseudo hydrostatic pressure [Pa]
    - p     = Static pressure [Pa]
    - h     = Height in the opposite direction to gravity
    - hRef  = Reference height in the opposite direction to gravity
    - ρ     = Density
    - g     = Acceleration due to gravity [m/s²]

    Parameters
    ----------
    value : float | int
        Value at the patch set uniformly
    p : float | int
        Total pressure
    rho : str
        Density field name (default `thermo:rho`)
    """
    def __init__(
            self,
            value: float,
            p: float,
            rho: str="thermo:rho"
        ):
        super().__init__(type="prghPressure", value=value)

        self.p = p
        self.rho = rho

    @property
    def p(self):
        return self._p

    @p.setter
    def p(self, p) -> None:
        check_type("p", p, (float, int))
        self._p = p
        self.__setitem__('p', f"uniform {p}")

    @property
    def rho(self):
        return self._rho

    @rho.setter
    def rho(self, rho) -> None:
        check_type("rho", rho, str)
        self._rho = rho
        self.__setitem__('rho', rho)
