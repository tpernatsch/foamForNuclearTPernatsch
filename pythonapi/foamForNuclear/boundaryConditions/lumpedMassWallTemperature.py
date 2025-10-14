from foamForNuclear.boundaryConditions.temperatureCoupledBase import TemperatureCoupledBase
from foamForNuclear.checkvalue import check_type


class LumpedMassWallTemperature(TemperatureCoupledBase):
    """
    Employs a lumped mass model for temperature.

    It considers a single temperature value for the whole patch and evaluates
    the temperature evolution using the net heat flux into the patch.

    Parameters
    ----------
    value : float | int
        Value at the patch set uniformly
    Cp : float | int
        Heat capacity [J/kg.K]
    mass : float | int
        Total mass [kg]
    kappaMethod : str
        Thermal conductivity option (default `fluidThermo`)
        The thermal conductivity kappa may be obtained by the following methods:

        - `lookup` = lookup volScalarField (or volSymmTensorField) with name \
            defined by `kappa`.
        - `fluidThermo` = use fluidThermo and default \
            `compressible::turbulenceModel` to calculate kappa.
        - `solidThermo` = use solidThermo kappa().
        - `directionalSolidThermo` = uses look up for volSymmTensorField for \
            transformed kappa vector. Field name definable in `alphaAni`, named \
            `Anialpha` in solid solver by default
        - `function` = kappa, alpha directly specified as Function1
        - `phaseSystem` = used for multiphase thermos

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
            value: float,
            Cp: float,
            mass: float,
            kappaMethod: str="fluidThermo",
            kappa: str=None,
            alpha: str=None,
            alphaAni: str=None,
            kappaValue=None,
            alphaValue=None
        ):
        super().__init__(
            type="lumpedMassWallTemperature",
            value=value,
            kappaMethod=kappaMethod,
            kappa=kappa,
            alpha=alpha,
            alphaAni=alphaAni,
            kappaValue=kappaValue,
            alphaValue=alphaValue
        )

        self.Cp = Cp
        self.mass = mass

    @property
    def Cp(self):
        return self._Cp

    @Cp.setter
    def Cp(self, Cp) -> None:
        check_type("Cp", Cp, (float, int))
        self._Cp = Cp
        self.__setitem__('Cp', Cp)

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, mass) -> None:
        check_type("mass", mass, (float, int))
        self._mass = mass
        self.__setitem__('mass', mass)
