from .boundaryCondition import Patch
from .fixedValue import FixedValue
from foamForNuclear.checkvalue import check_type


class TotalPressure(FixedValue):
    """
    The condition sets the static pressure at the patch pp based on a
    specification of the total pressure, p0. The mode of operation is determined
    via the input entries and the dimensions of the convective flux, phi (ϕ).

    Parameters
    ----------
    value : float | int
        Value at the patch set uniformly
    p0 : float | int
        Total pressure
    phi : str
        Flux field name (default: `phi`)
    U : str
        Velocity field name (default: `U`)
    rho : str
        Density field name  (default: `rho`)
    psi : str
        Compressibility field name (default: `none`)
    gamma : float
        Cp/Cv (default: `1`)
    """
    def __init__(
            self,
            p0: float,
            value: float=None,
            phi: str='phi',
            U: str='U',
            rho: str='rho',
            psi: str='none',
            gamma: float=1
        ):
        super().__init__(value=value)

        self.type = "totalPressure"
        self.phi = phi
        self.p0 = p0
        self.U = U
        self.rho = rho
        self.psi = psi
        self.gamma = gamma


    @property
    def p0(self):
        return self._p0

    @p0.setter
    def p0(self, p0) -> None:
        check_type("p0", p0, (float, int))
        self._p0 = p0
        self.__setitem__('p0', f"uniform {p0}")
        if (self.value is None):
            self.value = p0


    @property
    def phi(self):
        return self._phi

    @phi.setter
    def phi(self, phi) -> None:
        check_type("phi", phi, str)
        self._phi = phi
        self.__setitem__('phi', phi)


    @property
    def U(self):
        return self._U

    @U.setter
    def U(self, U) -> None:
        check_type("U", U, str)
        self._U = U
        self.__setitem__('U', U)


    @property
    def rho(self):
        return self._rho

    @rho.setter
    def rho(self, rho) -> None:
        check_type("rho", rho, str)
        self._rho = rho
        self.__setitem__('rho', rho)


    @property
    def psi(self):
        return self._psi

    @psi.setter
    def psi(self, psi) -> None:
        check_type("psi", psi, str)
        self._psi = psi
        self.__setitem__('psi', psi)


    @property
    def gamma(self):
        return self._gamma

    @gamma.setter
    def gamma(self, gamma) -> None:
        check_type("gamma", gamma, (float, int))
        self._gamma = gamma
        self.__setitem__('gamma', gamma)
