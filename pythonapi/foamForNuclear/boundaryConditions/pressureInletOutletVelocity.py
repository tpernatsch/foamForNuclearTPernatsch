from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type


class PressureInletOutletVelocity(Patch):
    """
    The `pressureInletOutletVelocity` is a velocity inlet/outlet boundary
    condition that applies a zero-gradient condition for outflow (as defined by
    the flux); and obtains velocity from the flux with specified inlet
    direction.

    The tangential patch velocity can be optionally specified.

    Parameters
    ----------
    value : float | int | Vector
        Value at the patch set uniformly
    phi : str (optional)
        Name of the flux field (default: `phi`)
    """

    def __init__(
            self,
            value: float,
            phi: str='phi'
        ):
        super().__init__(type="pressureInletOutletVelocity", value=value)

        self.phi = phi


    @property
    def phi(self):
        return self._phi

    @phi.setter
    def phi(self, phi) -> None:
        check_type("phi", phi, str, none_ok=True)
        self._phi = phi
        self.__setitem__('phi', phi)
