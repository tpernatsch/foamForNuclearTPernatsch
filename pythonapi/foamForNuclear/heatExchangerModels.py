from foamForNuclear.checkvalue import check_positive, check_type
from foamForNuclear.common import *


class HeatExchangerModel:
    """
    Model of heat exchanger to be used in phaseProperties

    Parameters
    ----------
    name : str
        Name of the heat exchanger
    primary : str
        Name of the cellZone acting as the primary circuit
    secondary: str
        Name of the cellZone acting as the secondary circuit
    volumetricArea: float
        Volumetric area in m2/m3
    wallConductance: float
        Wall conductance between the primary and secondary circuit
    """
    def __init__(
            self,
            name: str,
            primary: str,
            secondary: str,
            volumetricArea: float,
            wallConductance: float,
        ):
        self.name = name
        self.primary = primary
        self.secondary = secondary
        self.volumetricArea = volumetricArea
        self.wallConductance = wallConductance


    def __repr__(self, depth = 0):
        # Declare twice the Heat Exchanger, one per cellZone. This allows to
        # link 2 heat exchanegr side from 2 different mesh (e.g. 1 with sodium
        # fluid and 1 with water fluid)
        pDict = OpenFOAMDict({
            "primary": self.primary,
            "secondary": self.secondary,
            "volumetricArea": self.volumetricArea,
            "wallConductance": self.wallConductance
        })
        sDict = OpenFOAMDict({
            "primary": self.secondary,
            "secondary": self.primary,
            "volumetricArea": self.volumetricArea,
            "wallConductance": self.wallConductance
        })

        text = ""
        text += f"p{self.name}{pDict.__repr__(depth=depth)}"
        text += f"{depth*tab}s{self.name}{sDict.__repr__(depth=depth)}"
        return(text)


    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name) -> None:
        check_type("name", name, str)
        self._name = name


    @property
    def primary(self):
        return self._primary

    @primary.setter
    def primary(self, primary) -> None:
        check_type("primary", primary, str)
        self._primary = primary


    @property
    def secondary(self):
        return self._secondary

    @secondary.setter
    def secondary(self, secondary) -> None:
        check_type("secondary", secondary, str)
        self._secondary = secondary


    @property
    def volumetricArea(self):
        return self._volumetricArea

    @volumetricArea.setter
    def volumetricArea(self, volumetricArea) -> None:
        check_type("volumetricArea", volumetricArea, (int, float))
        check_positive("volumetricArea", volumetricArea)
        self._volumetricArea = volumetricArea


    @property
    def wallConductance(self):
        return self._wallConductance

    @wallConductance.setter
    def wallConductance(self, wallConductance) -> None:
        check_type("wallConductance", wallConductance, (int, float))
        self._wallConductance = wallConductance
