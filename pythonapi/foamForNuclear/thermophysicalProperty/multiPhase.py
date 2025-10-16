from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import List, OpenFOAMDict, addParameter
from foamForNuclear.openfoamFile import OpenFOAMFile
from .thermophysicalProperty import BaseThermophysicalProperty


class MultiPhase(BaseThermophysicalProperty):
    """
    Multi-phase model.

    Parameters
    ----------
    region : str
        Name of the region.
    phases : list[BaseThermophysicalProperty]
        List of phase as thermophysical property objects
    pMin : float
        Minimum pressure (default: 1e5 Pa)
    sigma : OpenFOAMDict
        Sigma model as an OpenFOAMDict
    """
    def __init__(
            self,
            region: str="",
            phases: list=[],
            pMin: float=1e5,
            sigma: OpenFOAMDict=OpenFOAMDict({
                'type': 'constant',
                'sigma': 0.07
            })
        ):
        super().__init__(region)

        self.phases: list[BaseThermophysicalProperty] = phases

        self.pMin = pMin

        self.sigma = sigma

        # self.transportProperties: TransportProperties = TransportProperties()


    @property
    def phases(self):
        return self._phases

    @phases.setter
    def phases(self, phases) -> None:
        check_type("phases", phases, list)
        self._phases = phases


    @property
    def pMin(self):
        return self._pMin

    @pMin.setter
    def pMin(self, pMin) -> None:
        check_type("pMin", pMin, (float, int))
        self._pMin = pMin


    @property
    def sigma(self):
        return self._sigma

    @sigma.setter
    def sigma(self, sigma) -> None:
        check_type("sigma", sigma, OpenFOAMDict)
        self._sigma = sigma


    def appendPhase(self, phase):
        check_type("phase", phase, BaseThermophysicalProperty)
        self.phases.append(phase)


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        # self.transportProperties.region = self.region

        # self.transportProperties.export_to_openfoam()

        phaseNames = []

        for phase in self.phases:
            phase.region = self.region
            phase.export_to_openfoam()
            phaseNames.append(phase.ext)

        text = ""

        text += addParameter("phases", List(phaseNames), isAddExtraLine=True)
        text += addParameter("pMin", self.pMin, isAddExtraLine=True)
        text += f"sigma{self.sigma}"

        return text
