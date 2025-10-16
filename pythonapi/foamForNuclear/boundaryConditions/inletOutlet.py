from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import Vector


class InletOutlet(Patch):
    """
    This boundary condition provides a generic outflow condition, with specified
    inflow for the case of return flow.

    Note: Sign conventions
    ----------------------
        - Positive flux (out of domain): apply zero-gradient condition
        - Negative flux (into of domain): apply the `inletValue` fixed-value

    Parameters
    ----------
    value : float | int
        Value at the patch set uniformly
    inletValue : float | int
        Inlet value for reverse flow
    phi : str
        Flux field name (default `phi`)
    """

    def __init__(
            self,
            value: float | Vector,
            inletValue: float | Vector,
            phi: str='phi'
        ):
        super().__init__(type="inletOutlet", value=value)

        self.phi = phi
        self.inletValue = inletValue


    @property
    def inletValue(self):
        return self._inletValue

    @inletValue.setter
    def inletValue(self, inletValue) -> None:
        check_type("inletValue", inletValue, (float, int, Vector))
        self._inletValue = inletValue
        self.__setitem__('inletValue', f"uniform {inletValue}")


    @property
    def phi(self):
        return self._phi

    @phi.setter
    def phi(self, phi) -> None:
        check_type("phi", phi, str)
        self._phi = phi
        self.__setitem__('phi', phi)
