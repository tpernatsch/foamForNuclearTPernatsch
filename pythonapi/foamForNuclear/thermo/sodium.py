import attr
from attr import define, field

from foamForNuclear.common import OpenFOAMDict, Polynome
from .thermophysicalProperty import BaseThermophysicalProperty, ThermoType

from foamForNuclear._attrs_tools import auto_type_validator


@define(
    slots=True,
    on_setattr=[
        attr.setters.convert,
        attr.setters.validate,
    ],
    field_transformer=auto_type_validator,
    repr=False,
)
class SodiumConst(BaseThermophysicalProperty):
    """
    Liquid Sodium with constant properties.
    """
    region: str = ""
    ext: str = ""

    pRef: float = 1e5
    molWeight: float = 22.989769
    rho: float = 1000
    Cp: float = 850
    Hf: float = 0
    Sf: float = 0
    mu: float = 1.8e-4
    kappa: float = 65
    Pr: float = 5.6664993e-3

    thermoType: ThermoType = field(init=False)

    def __attrs_post_init__(self):
        BaseThermophysicalProperty.__init__(self, self.region, ext=self.ext)

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


@define(
    slots=True,
    on_setattr=[
        attr.setters.convert,
        attr.setters.validate,
    ],
    field_transformer=auto_type_validator,
    repr=False,
)
class SodiumBoussinesq(BaseThermophysicalProperty):
    """
    Liquid Sodium using Boussinesq formalism.
    Default properties are those of liquid Sodium evaluated at 600 K.
    """
    region: str = ""
    ext: str = ""

    description: str = "Sodium properties"

    pRef: float | int = 1e5
    molWeight: float | int = 22.989769
    rho0: float | int = 1000
    T0: float | int = 600
    beta: float | int = 2.602975e-4
    Cp: float | int = 1301
    Hf: float | int = 0
    mu: float | int = 3.21e-4
    Pr: float | int = 5.6664993e-3

    thermoType: ThermoType = field(init=False)

    def __attrs_post_init__(self):
        BaseThermophysicalProperty.__init__(self, self.region, ext=self.ext)

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


@define(
    slots=True,
    on_setattr=[
        attr.setters.convert,
        attr.setters.validate,
    ],
    field_transformer=auto_type_validator,
    repr=False,
)
class SodiumPolynomial(BaseThermophysicalProperty):
    """
    Thermophysical properties based on
    `GeN-Foam/Tutorials/featureCases/1D_HX/onePhase/constant/fluidRegion/thermophysicalProperties`.

    In this properties, Cp is kept constant.
    """
    region: str = ""
    ext: str = ""

    description: str = (
        "rho: Discrepancy with experimental values below 1% in the\n"
        "     400 K - 1950 K range, deteriorates slowly\n"
        "mu:  Discrepancy with experimental values below 1% in the\n"
        "     400 K - 1400 K, deteriorates quickly after that\n"
        "kappa: Provided by evaluators for the 371 K - 2503.7 K\n"
    )

    pRef: float = 1e5
    molWeight: float = 22.989769
    rho: Polynome = field(factory=lambda: Polynome(1006.81, -2.17E-01, -2.83E-06, -6.54E-09))
    Cp: Polynome = field(factory=lambda: Polynome(1250))
    Hf: float = 0
    Sf: float = 0
    mu: Polynome = field(
        factory=lambda: Polynome(
            4.932462E-03,
            -2.654334E-05,
            6.527520E-08,
            -8.752858E-11,
            6.621667E-14,
            -2.657514E-17,
            4.402963E-21,
        )
    )
    kappa: Polynome = field(factory=lambda: Polynome(124.67, -1.1381e-1, 5.5226e-5, -1.1842e-8))

    thermoType: ThermoType = field(init=False)

    def __attrs_post_init__(self):
        BaseThermophysicalProperty.__init__(self, self.region, ext=self.ext)

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


@define(
    slots=True,
    on_setattr=[
        attr.setters.convert,
        attr.setters.validate,
    ],
    field_transformer=auto_type_validator,
    repr=False,
)
class SodiumVapourPerfectGas(BaseThermophysicalProperty):
    """
    Thermophysical properties of sodium vapour at 1 bar.
    From `1D_HX/twoPhase`.
    """
    region: str = ""
    ext: str = ""

    description: str = ""

    pRef: float = 100000
    molWeight: float = 22.989769
    Cp: float = 904.141
    Hf: float = 4217528
    Sf: float = 0
    mu: float = 1.8e-7
    Pr: float = 0.005

    thermoType: ThermoType = field(init=False)

    def __attrs_post_init__(self):
        BaseThermophysicalProperty.__init__(self, self.region, ext=self.ext)

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