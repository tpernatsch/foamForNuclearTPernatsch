from foamForNuclear.common import OpenFOAMDict, Polynome
from .thermophysicalProperty import BaseThermophysicalProperty, ThermoType


class SodiumConst(BaseThermophysicalProperty):
    """
    Liquid Sodium with constant properties.

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
        self.molWeight = 22.989769
        self.rho = 1000
        self.Cp = 850
        self.Hf = 0
        self.Sf = 0
        self.mu = 1.8e-4
        self.kappa = 65
        self.Pr = 5.6664993e-3

        self.thermoType = ThermoType(
            type="heRhoThermo",
            mixture="pureMixture",
            transport="const",
            thermo="hConst",
            equationOfState="rhoConst",
            specie="specie",
            energy="sensibleEnthalpy",
        )

    def update_mixture(self):
        self.mixture = OpenFOAMDict({
            "specie": OpenFOAMDict({
                "molWeight": self.molWeight
            }),
            "equationOfState": OpenFOAMDict({
                "rho": self.rho
            }),
            "thermodynamics": OpenFOAMDict({
                "Cp": self.Cp,
                "Hf": self.Hf,
                "Sf": self.Sf
            }),
            "transport": OpenFOAMDict({
                "mu": self.mu,
                "kappa": self.kappa,
                "Pr": self.Pr
            }),
        })


class SodiumBoussinesq(BaseThermophysicalProperty):
    """
    Liquid Sodium using Boussinesq formalism.

    Parameters
    ----------
    region : str
        Name of the region.
    ext : str
        Extension at the end of the file, e.g `".liquid"` (default `""`).
    """
    def __init__(self, region: str="", ext: str=""):
        super().__init__(region, ext=ext)

        self.description = "Properties are those of liquid Sodium evaluated at 600 K"

        self.pRef = 1e5
        self.molWeight = 22.989769
        self.rho0 = 1000
        self.T0 = 600
        self.beta = 2.602975e-4
        self.Cp = 1301
        self.Hf = 0
        self.mu = 3.21e-4
        self.Pr = 5.6664993e-3

        self.thermoType = ThermoType(
            type="heRhoThermo",
            mixture="pureMixture",
            equationOfState="Boussinesq",
            thermo="hConst",
            transport="const",
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


class SodiumPolynomial(BaseThermophysicalProperty):
    """
    Thermophysical properties based on `GeN-Foam/Tutorials/featureCases/1D_HX/onePhase/constant/fluidRegion/thermophysicalProperties`

    In this properties, Cp is kept constant

    Parameters
    ----------
    region : str
        Name of the region.
    ext : str
        Extension at the end of the file, e.g `".liquid"` (default `""`).
    """

    def __init__(self, region: str="", ext: str=""):
        super().__init__(region, ext=ext)

        self.description  = "rho: Discrepancy with experimental values below 1% in the\n"
        self.description += "     400 K - 1950 K range, deteriorates slowly\n"

        self.description += "mu:  Discrepancy with experimental values below 1% in the\n"
        self.description += "     400 K - 1400 K, deteriorates quickly after that\n"

        self.description += "kappa: Provided by evaluators for the 371 K - 2503.7 K\n"


        self.pRef = 1e5
        self.molWeight = 22.989769
        self.rho = Polynome(1006.81, -2.17E-01, -2.83E-06, -6.54E-09)
        self.Cp = Polynome(1250)
        self.Hf = 0
        self.Sf = 0
        self.mu = Polynome(4.932462E-03, -2.654334E-05, 6.527520E-08, -8.752858E-11, 6.621667E-14, -2.657514E-17, 4.402963E-21)
        self.kappa = Polynome(124.67, -1.1381e-1, 5.5226e-5, -1.1842e-8)

        self.thermoType = ThermoType(
            type="heRhoThermo",
            mixture="pureMixture",
            transport="polynomial",
            thermo="hPolynomial",
            equationOfState="icoPolynomial",
            specie="specie",
            energy="sensibleEnthalpy",
        )


    def update_mixture(self):
        self.mixture = OpenFOAMDict({
            "specie": OpenFOAMDict({
                "molWeight": self.molWeight
            }),
            "equationOfState": OpenFOAMDict({
                "rhoCoeffs<8>": self.rho
            }),
            "thermodynamics": OpenFOAMDict({
                "CpCoeffs<8>": self.Cp,
                "Hf": self.Hf,
                "Sf": self.Sf
            }),
            "transport": OpenFOAMDict({
                "muCoeffs<8>": self.mu,
                "kappaCoeffs<8>": self.kappa
            }),
        })


class SodiumVapourPerfectGas(BaseThermophysicalProperty):
    """
    Thermophysical properties of sodium vapour at 1 bar.
    From `1D_HX/twoPhase`
    """

    def __init__(self, region="", ext=""):
        super().__init__(region, ext=ext)

        self.description  = ""

        self.pRef = 100000
        self.molWeight = 22.989769
        self.Cp = 904.141
        self.Hf = 4217528
        self.Sf = 0
        self.mu = 1.8e-7
        self.Pr = 0.005

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
                "Hf": self.Hf,
                "Sf": self.Sf
            }),
            "transport": OpenFOAMDict({
                "mu": self.mu,
                "Pr": self.Pr
            }),
        })
