from foamForNuclear.boundaryConditions.boundaryCondition import Patch


class FixedFluxPressure(Patch):
    """
    This boundary condition sets the pressure gradient to the provided value
    such that the flux on the boundary is that specified by the velocity
    boundary condition.

    Parameters
    ----------
    value : float | int
        Value at the patch set uniformly
    """

    def __init__(self, value: float):
        super().__init__(type="fixedFluxPressure", value=value)
