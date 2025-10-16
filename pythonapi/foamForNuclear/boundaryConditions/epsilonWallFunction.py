from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type, check_value


_BLENDING_TYPES = {"stepwise", "max", "binomial", "exponential"}


class EpsilonWallFunction(Patch):
    """
    The epsilonWallFunction boundary condition provides a wall constraint on the
    turbulent kinetic energy dissipation rate, i.e. epsilon, and the turbulent
    kinetic energy production contribution, i.e. G, for low- and high-Reynolds
    number turbulence models.

    Parameters
    ----------
    value : float | int | Vector
        Value at the patch set uniformly
    lowReCorrection : bool
        Flag: apply low-Re correction
    blending : str
        Viscous/inertial sublayer blending method
            - `stepwise`   : Stepwise switch (discontinuous)
            - `max`        : Maximum value switch (discontinuous)
            - `binomial`   : Binomial blending (smooth)
            - `exponential`: Exponential blending (smooth)
    n : float | int
        Binomial blending exponent
    """

    def __init__(
            self,
            value: float,
            lowReCorrection: bool=False,
            blending: str='stepwise',
            n: float=2.0
        ):
        super().__init__(type="epsilonWallFunction", value=value)

        self.lowReCorrection = lowReCorrection
        self.blending = blending
        self.n = n

    def __repr__(self, depth = 0):
        self.__setitem__('lowReCorrection', self.lowReCorrection)
        self.__setitem__('blending', self.blending)
        self.__setitem__('n', self.n)
        self.__setitem__('value', f"uniform {self.value}")

        return super().__repr__(depth)


    @property
    def lowReCorrection(self):
        return self._lowReCorrection

    @lowReCorrection.setter
    def lowReCorrection(self, lowReCorrection) -> None:
        check_type("lowReCorrection", lowReCorrection, bool)
        self._lowReCorrection = lowReCorrection


    @property
    def blending(self):
        return self._blending

    @blending.setter
    def blending(self, blending) -> None:
        check_type("blending", blending, str)
        check_value("blending", blending, _BLENDING_TYPES)
        self._blending = blending


    @property
    def n(self):
        return self._n

    @n.setter
    def n(self, n) -> None:
        check_type("n", n, (float, int))
        self._n = n
