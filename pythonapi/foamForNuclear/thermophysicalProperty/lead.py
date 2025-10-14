from foamForNuclear.common import OpenFOAMDict
from .thermophysicalProperty import BaseThermophysicalProperty, ThermoType


class LeadBoussinesq(BaseThermophysicalProperty):
    """
    Thermophysical properties of lead from OECD. Handbook on Lead-bismuth
    Eutectic Alloy and Lead Properties, Materials Compatibility,
    Thermal-hydraulics and Technologies. Technical report, OECD/NEA Nuclear
    Science Committee, 2015.
    Table 1.1 Basic characteristics of reactor coolants.
    Table 2.19.1 Summary of the recommended correlations for main thermophysical
    properties of liquid lead (Pb).

    Parameters
    ----------
    region : str
        Name of the region.
    ext : str
        Extension at the end of the file, e.g `".liquid"` (default `""`).
    """
    def __init__(self, region: str="", ext: str=""):
        super().__init__(region, ext=ext)

        self.description  = "OECD. Handbook on Lead-bismuth Eutectic Alloy and Lead Properties,\n"
        self.description += "  Materials Compatibility, Thermal-hydraulics and Technologies. Technical\n"
        self.description += "  report, OECD/NEA Nuclear Science Committee, 2015\n"
        self.description += "  Table 1.1: Basic characteristics of reactor coolants\n"
        self.description += "  Table 2.19.1: Summary of the recommended correlations for main\n"
        self.description += "  thermophysical properties of liquid lead (Pb)"

        self.pRef = 1e5
        self.molWeight = 207
        self.rho0 = 11441
        self.T0 = 0
        self.beta = 1.2795/11441
        self.Cp = 1.459504e+02
        self.Hf = 0
        self.mu = 2.037739e-03
        self.Pr = 0.017443

        self.thermoType = ThermoType(
            type="heRhoThermo",
            mixture="pureMixture",
            transport="const",
            thermo="hConst",
            equationOfState="Boussinesq",
            specie="specie",
            energy="sensibleEnthalpy",
        )

    def update_mixture(self):
        self.mixture = OpenFOAMDict({
            "specie": OpenFOAMDict({
                "molWeight": self.molWeight
            }),
            "equationOfState": OpenFOAMDict({
                "rho0": self.rho0,
                "T0": self.T0,
                "beta": self.beta
            }),
            "thermodynamics": OpenFOAMDict({
                "Cp": self.Cp,
                "Hf": self.Hf,
            }),
            "transport": OpenFOAMDict({
                "mu": self.mu,
                "Pr": self.Pr
            }),
        })
