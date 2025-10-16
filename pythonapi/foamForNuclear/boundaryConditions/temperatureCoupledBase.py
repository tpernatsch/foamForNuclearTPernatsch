from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type, check_value


_KAPPA_METHOD_TYPES = {
    "lookup", "fluidThermo", "solidThermo", "directionalSolidThermo",
    "function", "phaseSystem"
}


class TemperatureCoupledBase(Patch):
    """
    Common functions used in temperature coupled boundaries.

    Parameters
    ----------
    value : float | int
        Value at the patch set uniformly
    kappaMethod : str
        Thermal conductivity option (default `fluidThermo`)
        The thermal conductivity kappa may be obtained by the following methods:

        - `lookup` : lookup volScalarField (or volSymmTensorField) with name \
            defined by `kappa`
        - `fluidThermo` : use fluidThermo and default \
            `compressible::turbulenceModel` to calculate kappa
        - `solidThermo` : use solidThermo kappa()
        - `directionalSolidThermo`: uses look up for volSymmTensorField for \
            transformed kappa vector. Field name definable in `alphaAni`, named \
            `Anialpha` in solid solver by default
        - `function` : kappa, alpha directly specified as Function1
        - `phaseSystem` : used for multiphase thermos

    kappa : str
        Name of thermal conductivity field (default `None`)
    alpha : str
        Name of thermal diffusivity field (default `None`)
    alphaAni : str
        Name of non-isotropic alpha (default `None`)
    kappaValue
        Function1 supplying kappa (default `None`)
    alphaValue
        Function1 supplying alpha (default `None`)
    """
    def __init__(
            self,
            type: str,
            value: float,
            kappaMethod: str="fluidThermo",
            kappa: str=None,
            alpha: str=None,
            alphaAni: str=None,
            kappaValue=None,
            alphaValue=None
        ):
        super().__init__(type=type, value=value)

        self.kappaMethod = kappaMethod
        self.kappa = kappa
        self.alpha = alpha
        self.alphaAni = alphaAni
        self.kappaValue = kappaValue
        self.alphaValue = alphaValue

    @property
    def kappaMethod(self):
        return self._kappaMethod

    @kappaMethod.setter
    def kappaMethod(self, kappaMethod) -> None:
        check_type("kappaMethod", kappaMethod, str)
        check_value("kappaMethod", kappaMethod, _KAPPA_METHOD_TYPES)
        self._kappaMethod = kappaMethod
        self.__setitem__('kappaMethod', kappaMethod)

    @property
    def kappa(self):
        return self._kappa

    @kappa.setter
    def kappa(self, kappa) -> None:
        check_type("kappa", kappa, str, none_ok=True)
        self._kappa = kappa
        if (kappa is not None):
            self.__setitem__('kappa', kappa)

    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, alpha) -> None:
        check_type("alpha", alpha, str, none_ok=True)
        self._alpha = alpha
        if (alpha is not None):
            self.__setitem__('alpha', alpha)

    @property
    def alphaAni(self):
        return self._alphaAni

    @alphaAni.setter
    def alphaAni(self, alphaAni) -> None:
        check_type("alphaAni", alphaAni, str, none_ok=True)
        self._alphaAni = alphaAni
        if (alphaAni is not None):
            self.__setitem__('alphaAni', alphaAni)

    @property
    def kappaValue(self):
        return self._kappaValue

    @kappaValue.setter
    def kappaValue(self, kappaValue) -> None:
        self._kappaValue = kappaValue
        if (kappaValue is not None):
            self.__setitem__('kappaValue', f"{kappaValue}")

    @property
    def alphaValue(self):
        return self._alphaValue

    @alphaValue.setter
    def alphaValue(self, alphaValue) -> None:
        self._alphaValue = alphaValue
        if (alphaValue is not None):
            self.__setitem__('alphaValue', f"{alphaValue}")
