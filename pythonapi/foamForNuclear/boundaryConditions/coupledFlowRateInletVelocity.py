from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import Vector


class CoupledFlowRateInletVelocity(Patch):
    """
    Velocity inlet boundary condition either correcting the extrapolated
    velocity or creating a uniform velocity field normal to the patch adjusted
    to match the specified flow rate.

    For a mass-based flux:

    - the flow rate should be provided in kg/s
    - if `rho` is "none" the flow rate is in m3/s
    - otherwise `rho` should correspond to the name of the density field
    - if the density field cannot be found in the database, the user must \
        specify the inlet density using the `rhoInlet` entry

    For a volumetric-based flux:

    - the flow rate is in m3/s

    Notes:

    - `rhoInlet` is required for the case of a mass flow rate, where the \
        density field is not available at start-up
    - The value is positive into the domain (as an inlet)
    - May not work correctly for transonic inlets
    - Strange behaviour with potentialFoam since the U equation is not solved

    Usage
    -----

    Example of the boundary condition specification for a volumetric flow rate::

        CoupledFlowRateInletVelocity(
            value=Vector(0, 0, 0),
            volumetricFlowRate="mdot_in",
            vFlowInit=0.2,
            extrapolateProfile=True,
        )

    Example of the boundary condition specification for a mass flow rate::

        CoupledFlowRateInletVelocity(
            value=Vector(0, 0, 0),
            massFlowRate="mdot_in",
            vFlowInit=0.2,
            rho="rho",
            rhoInlet=1.0,
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
    rhoInlet : float
        Inlet density, only used with `massFlowRate`
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
            rhoInlet: float=None,
            extrapolateProfile: bool=False,
        ):
        super().__init__(type="coupledFlowRateInletVelocity", value=value)

        self.volumetricFlowRate = volumetricFlowRate
        self.vFlowInit = vFlowInit

        self.massFlowRate = massFlowRate
        self.mFlowInit = mFlowInit
        self.rho = rho
        self.rhoInlet = rhoInlet
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
            if (self.rhoInlet is not None):
                self.__setitem__("rhoInlet", self.rhoInlet)

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
    def rhoInlet(self):
        return self._rhoInlet

    @rhoInlet.setter
    def rhoInlet(self, rhoInlet) -> None:
        check_type("rhoInlet", rhoInlet, (int, float), none_ok=True)
        self._rhoInlet = rhoInlet


    @property
    def extrapolateProfile(self):
        return self._extrapolateProfile

    @extrapolateProfile.setter
    def extrapolateProfile(self, extrapolateProfile) -> None:
        check_type("extrapolateProfile", extrapolateProfile, bool, none_ok=True)
        self._extrapolateProfile = extrapolateProfile
        self.__setitem__('extrapolateProfile', self.extrapolateProfile)
