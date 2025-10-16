from foamForNuclear.boundaryConditions.boundaryCondition import Patch


class NutkWallFunction(Patch):
    """
    The nutkWallFunction boundary condition provides a wall constraint on the
    turbulent viscosity, i.e. nut, based on the turbulent kinetic energy, i.e.
    k, for low- and high-Reynolds number turbulence models.

    Parameters
    ----------
    value : float | int
        Value at the patch set uniformly
    """

    def __init__(self, value: float):
        super().__init__(type="nutkWallFunction", value=value)
