from __future__ import annotations
from typing import Optional, Dict

from foamForNuclear.common import FoamForNuclearDict, OpenFOAMDict
from foamForNuclear._attrs_tools import ffn_define, _to_List_float
from typing import ClassVar
from attrs import field, Factory, setters
from foamForNuclear.common import OpenFOAMListDict, Table, List, CheckedList
from foamForNuclear.porous_medium.heat_transfer import onb

from . import properties
from .behaviour import swelling, relocation, densification, failure, phase_transition
from . import laws, damage


def _sync_Tref_and_alpha_hook(inst, attr, value):
    """
    on_setattr hook that keeps inst.Tref and inst.thermalExpansion.Tref in sync.

    - When Tref is set: push it into the current thermalExpansion model.
    - When thermalExpansion is set: ensure its Tref matches inst.Tref.
    """
    if attr.name == "Tref":
        # value is the new Tref
        te = getattr(inst, "thermalExpansion", None)
        if isinstance(te, properties.thermal_expansion.ThermalExpansion):
            te.Tref = value

    elif attr.name == "thermalExpansion":
        te = value  # post-conversion value if we pipe after setters.convert
        if isinstance(te, properties.thermal_expansion.ThermalExpansion):
            # align its Tref to the instance's current Tref (if already set)
            tref = getattr(inst, "Tref", None)
            # print(te.Tref)
            # print(tref)
            if tref is not None:
                te.Tref = tref

    # always return the value to let attrs continue
    return value


@ffn_define
class Material(FoamForNuclearDict):
    TYPE: ClassVar[str] = "none"
    density: properties.density.Density | None = None
    conductivity: properties.conductivity.Conductivity | None = None
    heatCapacity: properties.heat_capacity.HeatCapacity | None = None
    emissivity: properties.emissivity.Emissivity | None = None
    YoungModulus: properties.young_modulus.YoungModulus | None = None
    PoissonRatio: properties.poisson_ratio.PoissonRatio | None = None
    thermalExpansion: properties.thermal_expansion.ThermalExpansion | None = field(
        default=None,
        on_setattr=setters.pipe(
            setters.convert,      # run converter on set
            _sync_Tref_and_alpha_hook,
        ),
    )
    Tref: float | int = field(
        default=293.0,
        metadata={"ffn_internal": True},
        on_setattr=setters.pipe(
            setters.convert,      # if decorator uses convert-on-set
            _sync_Tref_and_alpha_hook,
        ),
    )
    # Damage model
    damage: damage.Damage | None = None
    # Constituive law
    constitutiveLaw: laws.ConstitutiveLaw = field(factory=laws.Elasticity)

    # Force sync between Tref and thermal expansion
    def __attrs_post_init__(self) -> None:
        self.Tref = self.Tref
        super().__attrs_post_init__()


# ------------ Constant Material and its functions  ---------
def _to_const_density(value):
    if isinstance(value, properties.density.Density):
        return value
    if isinstance(value, (int, float)):
        return properties.density.Constant(value=value)
    raise TypeError(f"density must be a Density model or a number, got {type(value)}")


def _to_const_conductivity(value) -> properties.conductivity.Conductivity:
    if isinstance(value, properties.conductivity.Conductivity):
        return value
    if isinstance(value, (int, float)):
        return properties.conductivity.Constant(value=value)
    raise TypeError(f"conductivity must be a Conductivity model or a number, got {type(value)}")


def _to_const_heat_capacity(value):
    if isinstance(value, properties.heat_capacity.HeatCapacity):
        return value
    if isinstance(value, (int, float)):
        return properties.heat_capacity.Constant(value=value)
    raise TypeError(f"heatCapacity must be a HeatCapacity model or a number, got {type(value)}")


def _to_const_emissivity(value):
    if isinstance(value, properties.emissivity.Emissivity):
        return value
    if isinstance(value, (int, float)):
        return properties.emissivity.Constant(value=value)
    raise TypeError(f"emissivity must be an Emissivity model or a number, got {type(value)}")


def _to_const_E(value):
    if isinstance(value, properties.young_modulus.YoungModulus):
        return value
    if isinstance(value, (int, float)):
        return properties.young_modulus.Constant(value=value)
    raise TypeError(f"YoungModulus must be a YoungModulus model or a number, got {type(value)}")


def _to_const_nu(value):
    if isinstance(value, properties.poisson_ratio.PoissonRatio):
        return value
    if isinstance(value, (int, float)):
        return properties.poisson_ratio.Constant(value=value)
    raise TypeError(f"PoissonRatio must be a PoissonRatio model or a number, got {type(value)}")


def _to_const_alpha(value):
    if isinstance(value, properties.thermal_expansion.ThermalExpansion):
        return value
    if isinstance(value, (int, float)):
        # pick a default Tref; you can make this configurable later
        return properties.thermal_expansion.Constant(value=value, Tref=293)
    raise TypeError(
        f"thermalExpansion must be a ThermalExpansion model or a number, got {type(value)}"
    )


@ffn_define
class Constant(Material):
    TYPE: ClassVar[str] = "constant"

    # override fields to add converters (types stay the same)
    density: float | properties.density.Density = field(
        default=1000, converter=_to_const_density)
    conductivity: float | properties.conductivity.Conductivity = field(
        default=20, converter=_to_const_conductivity)
    heatCapacity: float | properties.heat_capacity.HeatCapacity = field(
        default=200, converter=_to_const_heat_capacity)
    emissivity: float | properties.emissivity.Emissivity = field(
        default=0.8, converter=_to_const_emissivity)
    YoungModulus: float | properties.young_modulus.YoungModulus = field(
        default=100e9, converter=_to_const_E)
    PoissonRatio: float | properties.poisson_ratio.PoissonRatio = field(
        default=0.3, converter=_to_const_nu)
    thermalExpansion: float | properties.thermal_expansion.ThermalExpansion = field(
        default=1e-5,
        converter=_to_const_alpha,
        on_setattr=setters.pipe(
            setters.convert,      # run converter on set
            _sync_Tref_and_alpha_hook,
        ),
    )



# @ffn_define
# class Buffer(Material):
#     TYPE: ClassVar[str] = "buffer"
#     swelling: swelling.Swelling | None

    # theoreticalDensity: float | int = 2250.0

#     @classmethod
#     def preset(
#         cls,
#         *,
#         name: str,
#         initialDensity: float | int = 1000, theoreticalDensity: float | int = 2250,
#         initialConductivity: float | int = 0.5, theoreticalConductivity: float | int = 4,
#         cp: float | int =720, emissivity: float | int = 0.0,
#         nu: float | int = 0.33, Tref: float | int = 293
#     ) -> Constant:
#         return cls(
#             name=name,
#             density=properties.density.Constant(rho=initialDensity),
#             conductivity=properties.conductivity.BufferParfume(
#                 initialDensity=initialDensity, theoreticalDensity=theoreticalDensity,
#                 initialConductivity=initialConductivity, theoreticalConductivity=theoreticalConductivity
#             ),
#             heatCapacity=properties.heat_capacity.Constant(cp=cp),
#             emissivity=properties.emissivity.Constant(emissivity=emissivity),
#             YoungModulus=properties.young_modulus.BufferParfume(),
#             PoissonRatio=properties.poisson_ratio.Constant(nu=nu),
#             thermalExpansion=properties.thermal_expansion.BufferParfume(Tref=Tref),
#             swelling=swelling.BufferParfume()
#         )


# @ffn_define
# class Hastelloy(Material):
#     TYPE: ClassVar[str] = "hastelloy"
#     # density: properties.density.Density = field(factory=properties.density.Constant)
#     # conductivity: properties.conductivity.Conductivity = field(factory=properties.conductivity.HastelloyNSwindeman)
#     # heatCapacity: properties.heat_capacity.HeatCapacity = field(factory=properties.heat_capacity.Constant)
#     # emissivity: properties.emissivity.Emissivity = field(factory=properties.emissivity.Constant)
#     # YoungModulus: properties.young_modulus.YoungModulus = field(factory=properties.young_modulus.WatrousHastelloyN)
#     # PoissonRatio: properties.poisson_ratio.PoissonRatio = field(factory=properties.poisson_ratio.Constant)
#     # thermalExpansion: properties.thermal_expansion.ThermalExpansion = field(factory=properties.thermal_expansion.HastelloyNSwindeman)
#     # swelling: swelling.Swelling = field(factory=swelling.HastelloyNWrightSham)


# @ffn_define
# class Inconel600(Material):
#     TYPE: ClassVar[str] = "inconel600"
#     density: properties.density.Density = field(factory=properties.density.Constant)


# @ffn_define
# class Molybdenum(Material):
#     TYPE: ClassVar[str] = "molybdenum"
#     # density: properties.density.Density = field(factory=properties.density.MoConstant)
#     # conductivity: properties.conductivity.Conductivity = field(factory=properties.conductivity.Molybdenum)
#     # heatCapacity: properties.heat_capacity.HeatCapacity = field(factory=properties.heat_capacity.Molybdenum)
#     # emissivity: properties.emissivity.Emissivity = field(factory=properties.emissivity.MoConstant)
#     # YoungModulus: properties.young_modulus.YoungModulus = field(factory=properties.young_modulus.Molybdenum)
#     # PoissonRatio: properties.poisson_ratio.PoissonRatio = field(factory=properties.poisson_ratio.MoConstant)
#     # thermalExpansion: properties.thermal_expansion.ThermalExpansion = field(factory=properties.thermal_expansion.Molybdenum)
#     # swelling: swelling.Swelling = field(factory=swelling.FeCrAl)


# @ffn_define
# class PyC(Material):
#     TYPE: ClassVar[str] = "PyC"
#     # density: properties.density.Density = field(factory=properties.density.Constant)
#     # conductivity: properties.conductivity.Conductivity = field(factory=properties.conductivity.Constant)
#     # heatCapacity: properties.heat_capacity.HeatCapacity = field(factory=properties.heat_capacity.Constant)
#     # emissivity: properties.emissivity.Emissivity = field(factory=properties.emissivity.Emissivity)
#     # YoungModulus: properties.young_modulus.YoungModulus = field(factory=properties.young_modulus.PyCParfume)
#     # PoissonRatio: properties.poisson_ratio.PoissonRatio = field(factory=properties.poisson_ratio.Constant)
#     # thermalExpansion: properties.thermal_expansion.ThermalExpansion = field(factory=properties.thermal_expansion.PyCParfume)
#     # swelling: swelling.Swelling = field(factory=swelling.PyCParfume)


# @ffn_define
# class SiC(Material):
#     TYPE: ClassVar[str] = "SiC"
#     # density: properties.density.Density = field(factory=properties.density.Constant)
#     # conductivity: properties.conductivity.Conductivity = field(factory=properties.conductivity.SiCParfume)
#     # heatCapacity: properties.heat_capacity.HeatCapacity = field(factory=properties.heat_capacity.SiCSnead)
#     # emissivity: properties.emissivity.Emissivity = field(factory=properties.emissivity.Emissivity)
#     # YoungModulus: properties.young_modulus.YoungModulus = field(factory=properties.young_modulus.SiCParfume)
#     # PoissonRatio: properties.poisson_ratio.PoissonRatio = field(factory=properties.poisson_ratio.Constant)
#     # thermalExpansion: properties.thermal_expansion.ThermalExpansion = field(factory=properties.thermal_expansion.SiCParfume)


# @ffn_define
# class Steel1515Ti(Material):
#     TYPE: ClassVar[str] = "Steel1515Ti"
#     # density: properties.density.Density = field(factory=properties.density.Steel1515TiSchumann)
#     # conductivity: properties.conductivity.Conductivity = field(factory=properties.conductivity.Steel1515TiTobbe)
#     # heatCapacity: properties.heat_capacity.HeatCapacity = field(factory=properties.heat_capacity.Steel1515TiBanerjee)
#     # emissivity: properties.emissivity.Emissivity = field(factory=properties.emissivity.ZircaloyConstant)
#     # YoungModulus: properties.young_modulus.YoungModulus = field(factory=properties.young_modulus.Steel1515TiTobbe)
#     # PoissonRatio: properties.poisson_ratio.PoissonRatio = field(factory=properties.poisson_ratio.Steel1515TiTobbe)
#     # thermalExpansion: properties.thermal_expansion.ThermalExpansion = field(factory=properties.thermal_expansion.Steel1515TiGehr)
#     # swelling: swelling.Swelling = field(factory=swelling.Steel1515TiGeneralized)


# ---- Fuel materials and related classes ----
@ffn_define
class ActinideDict(FoamForNuclearDict):
    massNumbers: list[int]
    weightFractions: list[float]
    ratioOverMetal: float


@ffn_define
class FuelMaterial(Material):
    TYPE: ClassVar[str] = "fuel"

    # behavior models (off by default)
    densification: densification.Densification | None = None
    swelling:      swelling.Swelling         | None = None
    relocation:    relocation.Relocation     | None = None
    failure:       failure.Failure           | None = None

    # composition/spec
    densityFraction: float | int
    theoreticalDensity: float | int
    oxygenMetalRatio: float | int
    rGrain: float
    dishFraction: float | int
    GdContent: float | int

    isotopes: dict[str, ActinideDict] = field(factory=dict,)


@ffn_define
class UO2(FuelMaterial):
    TYPE: ClassVar[str] = "UO2"

    # property models
    density: properties.density.Density = field(
        factory=properties.density.UO2Constant)
    conductivity: properties.conductivity.Conductivity = field(
        factory=properties.conductivity.UO2Matpro)
    heatCapacity: properties.heat_capacity.HeatCapacity = field(
        factory=properties.heat_capacity.UO2Matpro)
    emissivity: properties.emissivity.Emissivity = field(
        factory=properties.emissivity.UO2Relap)
    YoungModulus: properties.young_modulus.YoungModulus = field(
        factory=properties.young_modulus.UO2Matpro)
    PoissonRatio: properties.poisson_ratio.PoissonRatio = field(
        factory=properties.poisson_ratio.UO2Constant)
    thermalExpansion: properties.thermal_expansion.ThermalExpansion = field(
        factory=properties.thermal_expansion.UO2Relap)

    # composition/spec
    densityFraction: float | int = 0.95
    theoreticalDensity: float | int = 10960.0
    oxygenMetalRatio: float | int = 2.0
    rGrain: float = 10e-6
    dishFraction: float | int = 0.0
    GdContent: float | int = 0.0


# @ffn_define
# class UPuO2(Material):
#     TYPE: ClassVar[str] = "UPuO2"

    # densityFraction: float | int = 0.945
    # theoreticalDensity: float | int = 10430.0
#     # density: properties.density.Density = field(factory=properties.density.UPuO2Constant)
#     # conductivity: properties.conductivity.Conductivity = field(factory=properties.conductivity.MaUPuO2Magni)
#     # heatCapacity: properties.heat_capacity.HeatCapacity = field(factory=properties.heat_capacity.UPuO2Matpro)
#     # emissivity: properties.emissivity.Emissivity = field(factory=properties.emissivity.UO2Relap)
#     # YoungModulus: properties.young_modulus.YoungModulus = field(factory=properties.young_modulus.UPuO2Matpro)
#     # PoissonRatio: properties.poisson_ratio.PoissonRatio = field(factory=properties.poisson_ratio.UPuO2Constant)
#     # thermalExpansion: properties.thermal_expansion.ThermalExpansion = field(factory=properties.thermal_expansion.UPuO2Matpro)
#     # densification: densification.Densification = field(factory=densification.UO2Frapcon)
#     # swelling: swelling.Swelling = field(factory=swelling.UO2Frapcon)
#     # relocation: relocation.Relocation = field(factory=relocation.UO2Frapcon)
#     # failure: failure.Failure = field(factory=failure.Failure)
#     # # pore_velocity: pore_velocity.PoreVelocity = field(factory=pore_velocity.PoreVelocity)
#     # isotopes: OpenFOAMListDict = field(default=Factory(lambda self:
#     #     OpenFOAMListDict(IsotopesDict, "isotopes",), takes_self=True))


@ffn_define
class Zircaloy(Material):
    TYPE: ClassVar[str] = "Zircaloy"

    # property models
    density: properties.density.Density = field(
        factory=properties.density.ZircaloyIaea)
    conductivity: properties.conductivity.Conductivity = field(
        factory=properties.conductivity.ZircaloyRelap)
    heatCapacity: properties.heat_capacity.HeatCapacity = field(
        factory=properties.heat_capacity.ZircaloyIaea)
    emissivity: properties.emissivity.Emissivity = field(
        factory=properties.emissivity.ZircaloyConstant)
    YoungModulus: properties.young_modulus.YoungModulus = field(
        factory=properties.young_modulus.ZircaloyMatpro)
    PoissonRatio: properties.poisson_ratio.PoissonRatio = field(
        factory=properties.poisson_ratio.ZircaloyConstant)
    thermalExpansion: properties.thermal_expansion.ThermalExpansion = field(
        factory=properties.thermal_expansion.ZircaloyMatpro)

    # Behavior
    swelling: swelling.Swelling | None = None
    phase_transition: phase_transition.PhaseTransition | None = None
    failure: failure.Failure | None = None