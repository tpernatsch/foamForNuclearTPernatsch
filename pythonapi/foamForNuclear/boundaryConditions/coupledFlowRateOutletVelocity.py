from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import Vector


class CoupledFlowRateOutletVelocity(Patch):
    """
    Velocity outlet boundary condition either correcting the extrapolated
    velocity or creating a uniform velocity field normal to the patch adjusted
    to match the specified flow rate

    For a mass-based flux:

    - the flow rate should be provided in kg/s
    - if `rho` is "none" the flow rate is in m3/s
    - otherwise `rho` should correspond to the name of the density field
    - if the density field cannot be found in the database, the user must \
        specify the outlet density using the `rhoOutlet` entry

    For a volumetric-based flux:

    - the flow rate is in m3/s

    Notes:

    - `rhoOutlet` is required for the case of a mass flow rate, where the \
        density field is not available at start-up
    - The value is positive into the domain (as an outlet)
    - May not work correctly for transonic outlets
    - Strange behaviour with potentialFoam since the U equation is not solved

    Usage
    -----

    Example of the boundary condition specification for a volumetric flow rate::

        CoupledFlowRateOutletVelocity(
            value=Vector(0, 0, 0),
            volumetricFlowRate="mdot_in",
            vFlowInit=0.2,
            extrapolateProfile=True,
        )

    Example of the boundary condition specification for a mass flow rate::

        CoupledFlowRateOutletVelocity(
            value=Vector(0, 0, 0),
            massFlowRate="mdot_in",
            vFlowInit=0.2,
            rho="rho",
            rhoOutlet=1.0,
            extrapolateProfile=True,
        )

    Parameters
    ----------
    volumetricFlowRate : str
        FMI port name for volumetric flow rate
    vFlowInit : float
        Volumetric flow rate [m3/s], only used with `volumetricFlowRate`
    massFlowRate : str
        FMI port name for mass flow rate
    mFlowInit : float
        Mass flow rate [kg/s], only used with `massFlowRate`
    rho : str
        Density field name, only used with `massFlowRate` (default `rho`)
    rhoOutlet : float
        Outlet density, only used with `massFlowRate`
    extrapolateProfile : bool (optional)
        Extrapolate velocity profile (default `False`)
    """

    def __init__(
            self,
            value: Vector,
            volumetricFlowRate: str=None,
            vFlowInit: float=None,
            massFlowRate: str=None,
            mFlowInit: float=None,
            rho: str=None,
            rhoOutlet: float=None,
            extrapolateProfile: bool=False,
        ):
        super().__init__(type="coupledFlowRateOutletVelocity", value=value)

        self.volumetricFlowRate = volumetricFlowRate
        self.vFlowInit = vFlowInit

        self.massFlowRate = massFlowRate
        self.mFlowInit = mFlowInit
        self.rho = rho
        self.rhoOutlet = rhoOutlet
        self.extrapolateProfile = extrapolateProfile


    def __repr__(self, depth=0):
        if (self.volumetricFlowRate is not None):
            self.__setitem__("volumetricFlowRate", self.volumetricFlowRate)
            self.__setitem__("vFlowInit", self.vFlowInit)

        elif (self.massFlowRate is not None):
            self.__setitem__("massFlowRate", self.massFlowRate)
            self.__setitem__("mFlowInit", self.mFlowInit)
            if (self.rho is not None):
                self.__setitem__("rho", self.rho)
            if (self.rhoOutlet is not None):
                self.__setitem__("rhoOutlet", self.rhoOutlet)

        self.__setitem__("extrapolateProfile", self.extrapolateProfile)

        return super().__repr__(depth)


    @property
    def volumetricFlowRate(self):
        return self._volumetricFlowRate

    @volumetricFlowRate.setter
    def volumetricFlowRate(self, volumetricFlowRate) -> None:
        check_type("volumetricFlowRate", volumetricFlowRate, str, none_ok=True)
        self._volumetricFlowRate = volumetricFlowRate


    @property
    def vFlowInit(self):
        return self._vFlowInit

    @vFlowInit.setter
    def vFlowInit(self, vFlowInit) -> None:
        check_type("vFlowInit", vFlowInit, (int, float), none_ok=True)
        self._vFlowInit = vFlowInit


    @property
    def massFlowRate(self):
        return self._massFlowRate

    @massFlowRate.setter
    def massFlowRate(self, massFlowRate) -> None:
        check_type("massFlowRate", massFlowRate, str, none_ok=True)
        self._massFlowRate = massFlowRate


    @property
    def mFlowInit(self):
        return self._mFlowInit

    @mFlowInit.setter
    def mFlowInit(self, mFlowInit) -> None:
        check_type("mFlowInit", mFlowInit, (int, float), none_ok=True)
        self._mFlowInit = mFlowInit


    @property
    def rho(self):
        return self._rho

    @rho.setter
    def rho(self, rho) -> None:
        check_type("rho", rho, str, none_ok=True)
        self._rho = rho


    @property
    def rhoOutlet(self):
        return self._rhoOutlet

    @rhoOutlet.setter
    def rhoOutlet(self, rhoOutlet) -> None:
        check_type("rhoOutlet", rhoOutlet, (int, float), none_ok=True)
        self._rhoOutlet = rhoOutlet


    @property
    def extrapolateProfile(self):
        return self._extrapolateProfile

    @extrapolateProfile.setter
    def extrapolateProfile(self, extrapolateProfile) -> None:
        check_type("extrapolateProfile", extrapolateProfile, bool, none_ok=True)
        self._extrapolateProfile = extrapolateProfile
        self.__setitem__('extrapolateProfile', self.extrapolateProfile)
