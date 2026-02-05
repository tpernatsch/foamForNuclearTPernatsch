from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.common import Vector

class FixedDisplacementZeroShear(Patch):
    """
    Boundary condition that enforces a fixed displacement in the normal direction
    while applying a zero-shear condition in the tangential directions.

    This condition is typically used to constrain the normal motion of a boundary
    while allowing tangential sliding without shear resistance.

    Options
    -------
    value : list | Vector
        Displacement value prescribed uniformly on the patch.
        (required: True)
    """
    def __init__(self, value: list | Vector):
        super().__init__(type="fixedDisplacementZeroShear", value=value)
