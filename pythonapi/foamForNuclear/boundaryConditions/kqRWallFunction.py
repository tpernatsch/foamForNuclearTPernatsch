from foamForNuclear.boundaryConditions.boundaryCondition import Patch


class KqRWallFunction(Patch):
    """
    The kqRWallFunction boundary condition provides a simple wrapper around the
    zero-gradient condition, which can be used for the turbulent kinetic energy,
    i.e. k, square-root of turbulent kinetic energy, i.e. q (e.g. in qZeta
    turbulence model), and Reynolds stress tensor fields, i.e. R (e.g. in LRR
    turbulence model), for the cases of high Reynolds number flow using wall
    functions.

    Parameters
    ----------
    value : float | int
        Value at the patch set uniformly
    """

    def __init__(self, value: float):
        super().__init__(type="kqRWallFunction", value=value)
