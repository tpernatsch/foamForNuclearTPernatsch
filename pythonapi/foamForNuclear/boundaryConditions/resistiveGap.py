from foamForNuclear.checkvalue import check_type, check_value
from .boundaryCondition import Patch
from foamForNuclear.common import Table, List


_REGION_COUPLED_TYPES = {"regionCoupledOFFBEAT"}


class ResistiveGap(Patch):
    """
    Coupled boundary condition for modeling a **fixed thermal gap conductance**
    between two solid bodies.

    In this model, the gap conductance is prescribed directly by the user and
    remains constant in time. This makes the condition suitable for simplified
    studies or for cases where the gap behaviour is known a priori and does not
    need to be coupled to contact, gas, or radiation models.

    This boundary condition must be applied symmetrically on a *pair* of
    `regionCoupledOFFBEAT` patches coupled via AMI
    (Arbitrary Mesh Interface).

    Options
    -------
    value : scalar
        Initial value of the patch temperature field.
        (required: True)

    coupled : bool
        Enable coupling between the two opposing patches.
        (default: True; required: False)

    patchType : word
        Type of the underlying fvPatch. Must be set to
        `regionCoupledOFFBEAT`.
        (default: regionCoupledOFFBEAT; required: False)

    kappa : word
        Name of the thermal conductivity field.
        (default: k; required: False)

    alpha : scalar
        Fixed gap conductance value.
        (default: 5000; required: False)
    """

    def __init__(
            self,
            value: float,
            coupled: bool=True,
            patchType: str="regionCoupledOFFBEAT",
            kappa: str="k",
            alpha: float | int =5000
        ):
        super().__init__(type="resistiveGap", value=value)

        self.patchType = patchType
        self.kappa = kappa
        self.coupled = coupled
        self.alpha = alpha


    def __repr__(self, depth = 0):
        self.__setitem__('patchType', self.patchType)
        self.__setitem__('kappa', self.kappa)
        self.__setitem__('coupled', self.coupled)
        self.__setitem__('alpha', f"uniform {self.alpha}")

        return super().__repr__(depth)


    @property
    def patchType(self):
        return self._patchType

    @patchType.setter
    def patchType(self, patchType) -> None:
        check_type("patchType", patchType, str)
        check_value("patchType", patchType, _REGION_COUPLED_TYPES)
        self._patchType = patchType


    @property
    def kappa(self):
        return self._kappa

    @kappa.setter
    def kappa(self, kappa) -> None:
        check_type("kappa", kappa, str)
        self._kappa = kappa

    @property
    def coupled(self):
        return self._coupled

    @coupled.setter
    def coupled(self, coupled) -> None:
        check_type("coupled", coupled, bool)
        self._coupled = coupled

    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, alpha) -> None:
        check_type("alpha", alpha, (int, float), none_ok=True)
        self._alpha = alpha
        if (alpha is not None):
            self.__setitem__('alpha', f"uniform {alpha}")