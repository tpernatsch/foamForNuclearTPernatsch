from foamForNuclear.common import OpenFOAMDict, Polynome
from .thermophysicalProperty import BaseThermophysicalProperty, ThermoType


class WaterPerfectFluid(BaseThermophysicalProperty):
    """
    Water as a perfect fluid.

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
        self.molWeight = 18.0
        self.rho0 = 1027
        self.R = 3000
        self.Cp = 4195
        self.Hf = 0
        self.mu = 3.645e-4
        self.Pr = 2.289

        self.thermoType = ThermoType(
            type="heRhoThermo",
            mixture="pureMixture",
            transport="const",
            thermo="hConst",
            equationOfState="perfectFluid",
            specie="specie",
            energy="sensibleInternalEnergy",
        )

    def update_mixture(self):
        self.mixture = OpenFOAMDict({
            "specie": OpenFOAMDict({
                "molWeight": self.molWeight
            }),
            "equationOfState": OpenFOAMDict({
                "R": self.R,
                "rho0": self.rho0
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

        self.transportProperties = OpenFOAMDict({
            'transportModel': 'Newtonian',
            'nu': 1e-06,
            'rho': 1000
        })


class WaterPolynomial(BaseThermophysicalProperty):
    """
    Thermophysical properties of water at 147 bars. rho, Cp, mu and kappa are
    defined as polynomes.

    Parameters
    ----------
    region : str
        Name of the region.
    ext : str
        Extension at the end of the file, e.g `".liquid"` (default `""`).
    """
    def __init__(self, region: str="", ext: str=""):
        super().__init__(region, ext=ext)

        self.description  = "rho: Derived from water data AT SATURATION in the 0.01-350 Celsius range\n"
        self.description += "     from NIST:\n"
        self.description += "     https://www.nist.gov/system/files/documents/srd/NISTIR5078-Tab1.pdf\n"

        self.description += "Cp:  Derived from water data AT SATURATION in the 200-350 Celsius range\n"
        self.description += "     from NIST:\n"
        self.description += "     https://www.nist.gov/system/files/documents/srd/NISTIR5078-Tab1.pdf\n"

        self.description += "mu:  Derived from water data AT SATURATION in the 200-350 Celsius range\n"
        self.description += "     From Nuclear Systems I: Thermal Hydraulic Fundamentals (No, Kazimi)\n"

        self.description += "kappa: Derived from water data AT SATURATION in the 0.01-350 Celsius range\n"
        self.description += "     From Nuclear Systems I: Thermal Hydraulic Fundamentals (No, Kazimi)\n"


        self.pRef = 14739394.95
        self.molWeight = 18.01528
        self.rho = Polynome(-5.387953E+03, 2.325700E+01, -2.197000E-02)
        self.Cp = Polynome(37571600.18, -459408.5397, 2324.788772, -6.235974159, 0.009357494, -7.45183E-06, 2.46163E-09)
        self.Hf = 0
        self.Sf = 0
        self.mu = Polynome(7.920000E-03, -3.991780E-05, 6.845860E-08, -3.952670E-11)
        self.kappa = Polynome(8.25128, -0.04077, 7.54E-05, -4.82E-08)

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

class WaterConst(BaseThermophysicalProperty):
    """
    Constant thermophysical properties of water at 14739394.95 Pa. Based on the
    `WaterPolynomial` thermophysical properties object.

    Parameters
    ----------
    T : float
        Temperature at which the thermophysical properties are computed.
    region : str
        Name of the region.
    ext : str
        Extension at the end of the file, e.g `".liquid"` (default `""`).
    """
    def __init__(self, T: float, region: str="", ext: str=""):
        super().__init__(region, ext)

        waterPolynomial = WaterPolynomial()

        self.pRef = waterPolynomial.pRef
        self.molWeight = waterPolynomial.molWeight
        self.rho = waterPolynomial.rho.value(T)
        self.Cp = waterPolynomial.Cp.value(T)
        self.Hf = waterPolynomial.Hf
        self.Sf = waterPolynomial.Sf
        self.mu = waterPolynomial.mu.value(T)
        self.kappa = waterPolynomial.kappa.value(T)
        self.Pr = 2.289

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


class WaterVapourPerfectGas(BaseThermophysicalProperty):
    """
    Thermophysical properties of water vapour at 147 bars.
    From `1D_CHF/imposedPower`

    Parameters
    ----------
    region : str
        Name of the region.
    ext : str
        Extension at the end of the file, e.g `".liquid"` (default `""`).
    """
    def __init__(self, region: str="", ext: str=""):
        super().__init__(region, ext=ext)

        self.description  = ""

        self.pRef = 14739394.95
        self.molWeight = 18.01528
        self.Cp = 10000
        self.Hf = 0
        self.Sf = 0
        self.mu = 3e-2
        self.Pr = 2

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
