from foamForNuclear.common import OpenFOAMDict
from .thermophysicalProperty import BaseThermophysicalProperty, ThermoType


class AirPerfectGas(BaseThermophysicalProperty):
    """
    Thermophysical properties of Air using perfect gas state equation.

    Parameters
    ----------
    region : str
        Name of the region.
    ext : str
        Extension at the end of the file, e.g `".liquid"` (default `""`).
    """
    def __init__(self, region: str="", ext: str=""):
        super().__init__(region, ext=ext)

        self.pRef = 1e5
        self.molWeight = 28.9
        self.Cp = 1007
        self.Hf = 0
        self.mu = 1.84e-05
        self.Pr = 0.7

        self.thermoType = ThermoType(
            type="heRhoThermo",
            mixture="pureMixture",
            transport="const",
            thermo="hConst",
            equationOfState="perfectGas",
            specie="specie",
            energy="sensibleInternalEnergy",
        )

        self.transportProperties = OpenFOAMDict({
            'transportModel': 'Newtonian',
            'nu': 1.48e-5,
            'rho': 1
        })

    def update_mixture(self):
        self.mixture = OpenFOAMDict({
            "specie": OpenFOAMDict({
                "molWeight": self.molWeight
            }),
            "thermodynamics": OpenFOAMDict({
                "Cp": self.Cp,
                "Hf": self.Hf
            }),
            "transport": OpenFOAMDict({
                "mu": self.mu,
                "Pr": self.Pr
            }),
        })
