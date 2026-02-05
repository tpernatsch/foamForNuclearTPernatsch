from foamForNuclear.common import OpenFOAMDict, Polynome
from .thermophysicalProperty import BaseThermophysicalProperty, ThermoType


class HydrogenPerfectGas(BaseThermophysicalProperty):
    """
    Hydrogen perfect gas model.

    Parameters
    ----------
    region : str
        Name of the region.
    ext : str
        Extension at the end of the file, e.g `".liquid"` (default `""`).
    """
    def __init__(self, region="", ext=""):
        super().__init__(region, ext=ext)

        self.description  = "Hydrogen values at 600K and 3.7 MPa pressure from NIST.\n"
        self.description += "At this temperature and pressure, the hydrogen is supercritical.\n"
        self.description += r"https://webbook.nist.gov/cgi/fluid.cgi?P=3.7&TLow=80&THigh=1000&TInc=10&Digits=5&ID=C1333740&Action=Load&Type=IsoBar&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm&RefState=DEF"

        self.pRef = 3.7e6
        self.molWeight = 2.016 # kg/kmol
        self.Cp = 14566 # J/kg/K
        self.Cv = 10439 # J/kg/K
        self.Hf = 0
        self.Sf = 0
        self.Tref = 600
        self.mu = 1.4462e-05 # Pa.s
        self.Pr = 0.701 # From https://lambdageeks.com/prandtl-number/

        self.thermoType = ThermoType(
            type="heRhoThermo",
            mixture="pureMixture",
            transport="const",
            thermo="hConst",
            equationOfState="perfectGas",
            specie="specie",
            energy="sensibleEnthalpy",
        )

    def update_mixture(self):
        self.mixture = OpenFOAMDict({
            "specie": OpenFOAMDict({
                "molWeight": self.molWeight
            }),
            "thermodynamics": OpenFOAMDict({
                "Cp": self.Cp,
                "Cv": self.Cv,
                "Hf": self.Hf,
                "Sf": self.Sf,
                "Tref": self.Tref,
            }),
            "transport": OpenFOAMDict({
                "mu": self.mu,
                "Pr": self.Pr
            }),
        })


class HydrogenPengRobinsonGas(BaseThermophysicalProperty):
    """
    Hydrogen thermophysical properties using the Peng-Robinson formalism.

    Parameters
    ----------
    region : str
        Name of the region.
    ext : str
        Extension at the end of the file, e.g `".liquid"` (default `""`).

    Attributes
    ----------
    Tc : float
        Critical temperature in K.
    Vc : float
        Critical volume in m3/kmol.
    Pc : float
        Critical pressure in Pa.
    omega : float
        Acentric factor.
    """

    def __init__(self, region: str="", ext: str=""):
        super().__init__(region, ext=ext)

        self.description  = "Hydrogen values at 600K and 3.7 MPa pressure from NIST.\n"
        self.description += "At this temperature and pressure, the hydrogen is supercritical.\n"
        self.description += r"https://webbook.nist.gov/cgi/fluid.cgi?P=3.7&TLow=80&THigh=1000&TInc=10&Digits=5&ID=C1333740&Action=Load&Type=IsoBar&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm&RefState=DEF"

        self.pRef = 1e5
        self.nMoles = 1
        self.molWeight = 2.016 # kg/kmol
        self.Tc = 33.145; # K, Critical temperature
        self.Vc = 1.0/(31.262/self.molWeight) # m3/kmol, Critical volume
        self.Pc = 1.2964e6 # Pa, Critical pressure
        self.omega = -0.219; # Acentric factor
        self.Cp = Polynome(
            1.266445630e+04, 1.17094171e+01, -2.98179826e-02, 3.81902598e-05,
            -2.54178633e-08, 9.28566400e-12, -1.77480852e-15, 1.39051574e-19
        )
        self.Hf = 0
        self.Sf = 0
        self.Tref = 600
        self.mu = Polynome(8.9385e-6)
        self.kappa = Polynome(0.01238494, 0.00057219)
        self.Pr = 0.701 # From https://lambdageeks.com/prandtl-number/

        self.thermoType = ThermoType(
            type="heRhoThermo",
            mixture="pureMixture",
            transport="polynomial",
            thermo="hPolynomial",
            equationOfState="PengRobinsonGas",
            specie="specie",
            energy="sensibleEnthalpy",
        )

    def update_mixture(self):
        self.mixture = OpenFOAMDict({
            "specie": OpenFOAMDict({
                "nMoles": self.nMoles,
                "molWeight": self.molWeight
            }),
            "equationOfState": OpenFOAMDict({
                "Tc": self.Tc,
                "Vc": self.Vc,
                "Pc": self.Pc,
                "omega": self.omega
            }),
            "thermodynamics": OpenFOAMDict({
                "CpCoeffs<8>": self.Cp,
                "Hf": self.Hf,
                "Sf": self.Sf,
                "Tref": self.Tref,
            }),
            "transport": OpenFOAMDict({
                "muCoeffs<8>": self.mu,
                "kappaCoeffs<8>": self.kappa,
                "Pr": self.Pr
            }),
        })


class Hydrogen(BaseThermophysicalProperty):
    """
    Hydrogen model that model the molecular hydrogen dissociation at
    high-temperature.

    Parameters
    ----------
    region : str
        Name of the region.
    ext : str
        Extension at the end of the file, e.g `".liquid"` (default `""`).
    """

    def __init__(self, region: str="", ext: str=""):
        super().__init__(region, ext=ext)

        self.description =  "Hydrogen model that model the molecular hydrogen dissociation at\n"
        self.description += "high-temperature"

        self.thermoType = ThermoType(
            type="heRhoThermo",
            mixture="pureMixture",
            properties="liquid",
            energy="sensibleEnthalpy",
        )

        self.mixture = OpenFOAMDict({
            "H2": None
        })