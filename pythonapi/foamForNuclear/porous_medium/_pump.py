from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import Vector
from ._phases import Structure
from foamForNuclear.timeProfile import TimeProfile


class Pump(Structure):
    """
    Allow defining a momentum source, including options for a time
    dependent source, and a source coming from an FMU.

    Parameters
    ----------
    zones : list[str]
        List of zone names to apply the momentum source
    volumeFraction : float
        Volume fraction of structure
    Dh : float
        Hydraulic diameter
    momentumSource: Vector
        Momentum source strength in kg.m/s / m3
    momentumSourceTimeProfile : TimeProfile, optional
        Momentum source scaling factor. Can be used with table or via external
        FMI input.
    """

    def __init__(
            self,
            zones: list[str]=[],
            volumeFraction: float=None,
            Dh: float=None,
            momentumSource: Vector=None,
            momentumSourceTimeProfile: TimeProfile=None
        ):
        super().__init__(zones, volumeFraction, Dh)

        self.momentumSource = momentumSource
        self.momentumSourceTimeProfile: TimeProfile = momentumSourceTimeProfile


    def __repr__(self, depth=0):
        if (self.momentumSourceTimeProfile is not None):
            self.__setitem__("momentumSourceTimeProfile", self.momentumSourceTimeProfile)

        return super().__repr__(depth)


    @property
    def momentumSource(self):
        return self._momentumSource

    @momentumSource.setter
    def momentumSource(self, momentumSource) -> None:
        check_type("momentumSource", momentumSource, Vector, none_ok=True)
        self._momentumSource = momentumSource
        self.__setitem__("momentumSource", momentumSource)


    @property
    def momentumSourceTimeProfile(self):
        return self._momentumSourceTimeProfile

    @momentumSourceTimeProfile.setter
    def momentumSourceTimeProfile(self, momentumSourceTimeProfile) -> None:
        check_type("momentumSourceTimeProfile", momentumSourceTimeProfile, TimeProfile, none_ok=True)
        self._momentumSourceTimeProfile = momentumSourceTimeProfile
