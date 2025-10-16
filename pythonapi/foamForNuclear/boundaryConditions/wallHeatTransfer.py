from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type


class WallHeatTransfer(Patch):
    """
    This boundary condition provides an enthalpy condition for wall heat
    transfer.

    Parameters
    ----------
    value : float | int
        Value at the patch set uniformly
    Tinf : float | int
        Wall temperature
    alphaWall : float | int
        Thermal diffusivity
    """
    def __init__(
            self,
            value: float,
            Tinf: float,
            alphaWall: float
        ):
        super().__init__(type="wallHeatTransfer", value=value)

        self.Tinf = Tinf
        self.alphaWall = alphaWall

    @property
    def Tinf(self):
        return self._Tinf

    @Tinf.setter
    def Tinf(self, Tinf) -> None:
        check_type("Tinf", Tinf, (float, int))
        self._Tinf = Tinf
        self.__setitem__('Tinf', f"uniform {Tinf}")

    @property
    def alphaWall(self):
        return self._alphaWall

    @alphaWall.setter
    def alphaWall(self, alphaWall) -> None:
        check_type("alphaWall", alphaWall, (float, int))
        self._alphaWall = alphaWall
        self.__setitem__('alphaWall', f"uniform {alphaWall}")
