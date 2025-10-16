from foamForNuclear.boundaryConditions.boundaryCondition import Patch


class MovingWallVelocity(Patch):
    """
    This boundary condition provides a velocity condition for cases with moving
    walls.

    Parameters
    ----------
    value : float | int | Vector
        Value at the patch set uniformly
    """

    def __init__(self, value: float):
        super().__init__(type="movingWallVelocity", value=value)
