from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type


class FlowRateInletVelocity(Patch):
    """
    Flow rate inlet velocity.

    Parameters
    ----------
    value : float | int
        Value at the patch set uniformly
    massFlowRate : float | int
        Mass flow rate value.
    volumetricFlowRate : float | int
        Volumetrc flow rate value.
    extrapolateProfile : bool
        A uniform plug flow is set by default. To set the profile according to
        the downstream cells, set `extrapolateProfile` to `True`
        (default `None`).
    rho : str
        Density field name (default `thermo:rho`)
    rhoInlet : float
        Inlet density
    """
    def __init__(
            self,
            value: float,
            massFlowRate: float=None,
            volumetricFlowRate: float=None,
            extrapolateProfile: bool=None,
            rho: str="thermo:rho",
            rhoInlet: float=None
        ):
        super().__init__(type="flowRateInletVelocity", value=value)

        self.massFlowRate = massFlowRate
        self.volumetricFlowRate = volumetricFlowRate
        self.extrapolateProfile = extrapolateProfile
        self.rho = rho
        self.rhoInlet = rhoInlet

    def __repr__(self, depth=0):
        if (
            (self.massFlowRate is None and self.volumetricFlowRate is None)
            or
            (self.massFlowRate is not None and self.volumetricFlowRate is not None)
        ):
            msg = "Please provide either 'massFlowRate' or 'volumetricFlowRate'"
            raise ValueError(msg)

        return super().__repr__(depth)

    @property
    def massFlowRate(self):
        return self._massFlowRate

    @massFlowRate.setter
    def massFlowRate(self, massFlowRate) -> None:
        check_type("massFlowRate", massFlowRate, (float, int), none_ok=True)
        self._massFlowRate = massFlowRate
        if (massFlowRate is not None):
            self.__setitem__('massFlowRate', massFlowRate)

    @property
    def volumetricFlowRate(self):
        return self._volumetricFlowRate

    @volumetricFlowRate.setter
    def volumetricFlowRate(self, volumetricFlowRate) -> None:
        check_type("volumetricFlowRate", volumetricFlowRate, (float, int), none_ok=True)
        self._volumetricFlowRate = volumetricFlowRate
        if (volumetricFlowRate is not None):
            self.__setitem__('volumetricFlowRate', volumetricFlowRate)

    @property
    def extrapolateProfile(self):
        return self._extrapolateProfile

    @extrapolateProfile.setter
    def extrapolateProfile(self, extrapolateProfile) -> None:
        check_type("extrapolateProfile", extrapolateProfile, bool, none_ok=True)
        self._extrapolateProfile = extrapolateProfile
        if (extrapolateProfile is not None):
            self.__setitem__('extrapolateProfile', extrapolateProfile)

    @property
    def rho(self):
        return self._rho

    @rho.setter
    def rho(self, rho) -> None:
        check_type("rho", rho, str)
        self._rho = rho
        self.__setitem__('rho', rho)

    @property
    def rhoInlet(self):
        return self._rhoInlet

    @rhoInlet.setter
    def rhoInlet(self, rhoInlet) -> None:
        check_type("rhoInlet", rhoInlet, (float, int), none_ok=True)
        self._rhoInlet = rhoInlet
        if (rhoInlet is not None):
            self.__setitem__('rhoInlet', rhoInlet)
