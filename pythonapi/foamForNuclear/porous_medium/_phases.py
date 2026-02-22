from types import NoneType
import numpy as np
from typing import Optional
from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import FoamForNuclearDict, List, OpenFOAMDict, Vector
from foamForNuclear.porous_medium.dispersed_diameter import DispersedDiameterModel
from foamForNuclear.porous_medium.power_models import PowerModel

from foamForNuclear.thermo import BaseThermophysicalProperty
from foamForNuclear.turbulence import TurbulenceProperties

from foamForNuclear._attrs_tools import ffn_define
from attrs import field, validators as v


_LATTICE_TYPES = {"square", "hexagon"}
_STATE_OF_MATTER_TYPES = {"liquid", "gas"}


@ffn_define
class PassiveProperties(FoamForNuclearDict):
    volumetricArea: int | float | None = None
    T: int | float | None = None
    rho: int | float | None = None
    Cp: int | float | None = None
    rhoCp: int | float | None = None


def _trigger_compute_hydraulic_parameters(instance, attribute, newValue):
    object.__setattr__(instance, attribute.name, newValue)
    print("trigger")
    instance.compute_hydraulic_parameters(
        pitch=instance.pitch,
        elementDiameter=instance.elementDiameter,
        latticeType=instance.latticeType,
        spacerThickness=instance.spacerThickness,
        gridFraction=instance.gridFraction,
        wireDiameter=instance.wireDiameter
    )


@ffn_define
class Structure(FoamForNuclearDict):
    """
    Structure property

    Two modes:
    - `volumeFraction` and `Dh` are provided, nothing happen
    - `pitch`, `elementDiameter`, `latticeType`, `spacerThickness`,
    `gridFraction` and `wireDiameter` are provided, `volumeFraction` and `Dh`
    are computed based on the previous parameters

    Parameters
    ----------
    zones : list[str]
        List of zone name to apply the structure model
    volumeFraction : float
        Volume fraction of structure
    Dh : float
        Hydraulic diameter
    pitch : str
        Center to center distance
    elementDiameter : float
        Element diameter
    latticeType : str
        Lattice type (`square`, `hexagon`)
    spacerThickness : float
        Spacer thickness (default 0)
    gridFraction : float
        Grid fraction defined as delta = (Total axial length of grid spacer)
        / (Axial length of the fuel bundle). Default 0
    wireDiameter : float
        Diameter of the wire (default 0) used in hexagonal assemblies
    localX : Vector
        Local X axes of the reference frame in this cellZone. These entries ARE
        NOT MANDATORY and default to `Vector(1, 0, 0)`.
    localZ : Vector
        Local Y axes of the reference frame in this cellZone. These entries ARE
        NOT MANDATORY and default to `Vector(0, 0, 1)`.
    localDhAnisotropy : Vector
        Used to compute a "directional hydraulic diameter". It builds upon
        the initial idea of GeN-Foam to use a directional Reyonlds number
        (i.e. different for reference frame direction). In particular, the
        i-th directional Reynolds will "see" a hydraulic diameter of
        Dh*localDhAnisotropy[i]. localDhAnisotropy is rotated to the local
        reference frame. This entry IS NOT MANDATORY and defaults to
        `Vector(1, 1, 1)`
    localTortuosity : Vector
        Tortuosity of the porous structure (in this cellZone) expressed in
        the local reference frame. This vector quantifies the deformation
        of diffusion lines (i.e. the path along which a quantity diffues)
        for each of the tree local reference frame axes. To provide an
        example, for a fluid in a bundle of pins, the tortuosity will be 1
        along the bundle axis and < 1 for directions transversal to the
        bundle axis. Components greater than 1 (not physical) mean that
        diffusion is easier along said direction. For components smaller
        than 1, the opposite holds.
        It is important to state that this ONLY IMPACTS THE DIFFUSION of
        fluid heat and momentum, not their advection. This entry
        IS NOT REQUIRED and deafults to `Vector(1, 1, 1)`

    Attributes
    ----------
    zones : list[str]
        List of zone name to apply the structure model
    volumeFraction : float
        Volume fraction of structure
    Dh : float
        Hydraulic diameter
    pitch : str
        Center to center distance
    elementDiameter : float
        Element diameter
    latticeType : str
        Lattice type (`square`, `hexagon`)
    spacerThickness : float
        Spacer thickness (default 0)
    gridFraction : float
        Grid fraction defined as delta = (Total axial length of grid spacer)
        / (Axial length of the fuel bundle). Default 0
    wireDiameter : float
        Diameter of the wire (default 0) used in hexagonal assemblies
    localX : Vector
        Local X axes of the reference frame in this cellZone. These entries ARE
        NOT MANDATORY and default to `Vector(1, 0, 0)`.
    localZ : Vector
        Local Y axes of the reference frame in this cellZone. These entries ARE
        NOT MANDATORY and default to `Vector(0, 0, 1)`.
    localDhAnisotropy : Vector
        Used to compute a "directional hydraulic diameter". It builds upon
        the initial idea of GeN-Foam to use a directional Reyonlds number
        (i.e. different for reference frame direction). In particular, the
        i-th directional Reynolds will "see" a hydraulic diameter of
        Dh*localDhAnisotropy[i]. localDhAnisotropy is rotated to the local
        reference frame. This entry IS NOT MANDATORY and defaults to
        `Vector(1, 1, 1)`
    localTortuosity : Vector
        Tortuosity of the porous structure (in this cellZone) expressed in
        the local reference frame. This vector quantifies the deformation
        of diffusion lines (i.e. the path along which a quantity diffues)
        for each of the tree local reference frame axes. To provide an
        example, for a fluid in a bundle of pins, the tortuosity will be 1
        along the bundle axis and < 1 for directions transversal to the
        bundle axis. Components greater than 1 (not physical) mean that
        diffusion is easier along said direction. For components smaller
        than 1, the opposite holds.
        It is important to state that this ONLY IMPACTS THE DIFFUSION of
        fluid heat and momentum, not their advection. This entry
        IS NOT REQUIRED and deafults to `Vector(1, 1, 1)`
    """
    zones: list | List = field(factory=list, metadata={"ffn_internal": True})
    volumeFraction: int | float | None = field(
        default=None,
        validator=v.optional(v.and_(v.instance_of((int, float)), v.ge(0), v.le(1))),
    )
    Dh: int | float | None = field(
        default=None,
        validator=v.optional(v.and_(v.instance_of((int, float)), v.ge(0))),
    )
    pitch: int | float | None = field(
        default=0,
        validator=v.optional(v.and_(v.instance_of((int, float)), v.ge(0),)),
        on_setattr=_trigger_compute_hydraulic_parameters
    )
    elementDiameter: int | float | None = field(
        default=0,
        validator=v.optional(v.and_(v.instance_of((int, float)), v.ge(0),)),
        on_setattr=_trigger_compute_hydraulic_parameters
    )
    latticeType: str | None = field(
        default="hexagon",
        validator=v.optional(v.and_(v.instance_of(str), v.in_(_LATTICE_TYPES),)),
        on_setattr=_trigger_compute_hydraulic_parameters
    )
    spacerThickness: int | float | None = field(
        default=0,
        validator=v.optional(v.and_(v.instance_of((int, float)), v.ge(0),)),
        on_setattr=_trigger_compute_hydraulic_parameters
    )
    gridFraction: int | float | None = field(
        default=0,
        validator=v.optional(v.and_(v.instance_of((int, float)), v.ge(0), v.le(1),)),
        on_setattr=_trigger_compute_hydraulic_parameters
    )
    wireDiameter: int | float | None = field(
        default=0,
        validator=v.optional(v.and_(v.instance_of((int, float)), v.ge(0),)),
        on_setattr=_trigger_compute_hydraulic_parameters
    )
    localX: Vector | None = None
    localZ: Vector | None = None
    localDhAnisotropy: Vector | None = None
    localTortuosity: Vector | None = None
    powerModel: PowerModel | None = None
    passiveProperties: PassiveProperties | None = None

    def __repr__(self, depth = 0):
        textZones = "\"" + ':'.join(self.zones) + "\""

        # if (self.powerModel is not None):
        #     self.__setitem__('powerModel', self.powerModel)

        # if (self.passiveProperties is not None):
        #     self.__setitem__('passiveProperties', self.passiveProperties)

        return textZones + super().__repr__(depth)

    def compute_hydraulic_parameters(
            self,
            pitch: float,
            elementDiameter: float,
            latticeType: str,
            spacerThickness: float=0,
            gridFraction: float=0,
            wireDiameter: float=0
        ):
        """
        Compute hydraulic diameter and volume fraction of structure based on
        N. Todreas, M. Kazimi and M. Massoud, "Nuclear Systems Volume 2,
        Elements of Thermal Hydraulic Design", Appendix J, page 603-606,
        CRV Press, 2022.

        Parameters
        ----------
        pitch : str
            Center to center distance
        elementDiameter : float
            Element diameter
        latticeType : str
            Lattice type (`square`, `hexagon`)
        spacerThickness : float
            Spacer thickness (default 0)
        gridFraction : float
            Grid fraction defined as delta = (Total axial length of grid spacer)
            / (Axial length of the fuel bundle). Default 0
        wireDiameter : float
            Diameter of the wire (default 0) used in hexagonal assemblies
        """
        print("hello")
        if (
            pitch is None
            or elementDiameter is None
            or latticeType is None
            or spacerThickness is None
            or gridFraction is None
            or wireDiameter is None
        ):
            return

        check_type("pitch", pitch, (float, int))
        check_type("elementDiameter", elementDiameter, (float, int))
        check_type("latticeType", latticeType, str)
        check_value("latticeType", latticeType, _LATTICE_TYPES)
        check_type("spacerThickness", spacerThickness, (float, int))
        check_type("gridFraction", gridFraction, (float, int))

        p = pitch
        D = elementDiameter
        t = spacerThickness
        delta = gridFraction
        Ds = wireDiameter

        if (latticeType == 'square'):
            # Total area
            At = p**2
            # Area for flow without spacer
            Af1 = At - np.pi*(D**2)/4
            # Average area for tlow with spacer
            Af1_s = Af1 - (2*p*t - t**2) * delta
            # Wetted perimeter for flow without spacer and duct
            Pw1 = np.pi * D
            # Average wetted perimeter including spacer but without duct
            Pw1_s = Pw1 + 4*(p-t)*delta

            self.Dh = 4*Af1_s / Pw1_s
            self.volumeFraction = 1 - Af1_s / At

        elif (latticeType == 'hexagon'):
            # Total area
            At = np.sqrt(3)/4 * p**2
            # Area for flow without wire wrap spacers
            Af1 = At - np.pi*(D**2) / 8
            # Area for flow including wire wrap spacers
            Af1_s = Af1 - np.pi*(Ds**2) / 8
            # Wetted perimeter including wire wrap spacers
            Pw1_s = 0.5 * np.pi * (D + Ds)

            self.Dh = 4*Af1_s / Pw1_s
            self.volumeFraction = 1 - Af1_s / At


@ffn_define
class Fluid(FoamForNuclearDict):
    """
    The presence of these dictionaries IS MANDATORY for two-phase simulations

    Parameters
    ----------
    stateOfMatter : str
        Supported entries are either `gas` or `liquid`; This entry is used by
        some drag, heat or mass transfer models as knowledge of which phase is
        gaseous(/vapourous) and which phase is liquid is necessary sometimes.
        This entry IS NOT MANDATORY in principle but might be REQUIRED by
        specific choices of models;
    thermoResidualAlpha : float
        Residual fluid volumeFraction below which the fluid temperature is not
        obtained from the enthalpy equation solution, but is set to the
        fluid1-fluid2 interfacial temperature. This is applied on a cell-by-cell
        basis (e.g. only in those cells whose volumeFraction of this fluid is
        below thermoResidualAlpha). This is meant to be used to stabilize the
        temperature field of a fluid being produced during phase change (either
        boiling or condensation). In fact, small inaccuracies in the calculation
        of the absolute enthalpy of a phase with a small volumeFraction will
        translate into large inaccuracies in the temperature. Note that
        thermoResidualAlpha is expressed relatively to the available volume for
        fluid flow, thus accounting for possible structures. This entry IS NOT
        MANDATORY and defaults to `0`;
    residualAlpha : float
        Residual fluid volumeFraction used to stabilize  equations when/if its
        volumeFraction tends to 0. This defaults to `1e-9` and IS NOT a mandatory
        entry;
    writeRestartFields : bool
        If true, write additional fields to disk that quickly allow residuals to
        converge to their pre-restart values if restarting the simulation from a
        certain time-step. This entry IS NOT MANDATORY and defaults to `False` to
        save disk space;
    """
    name: str = field(default="fluid", metadata={"ffn_internal": True})
    stateOfMatter: str = field(
        default="liquid",
        validator=v.optional(v.and_(v.instance_of(str), v.in_(_STATE_OF_MATTER_TYPES),)),
    )
    thermophysicalProperties: BaseThermophysicalProperty = field(factory=BaseThermophysicalProperty, metadata={"ffn_internal": True})
    turbulenceProperties: TurbulenceProperties = field(factory=TurbulenceProperties, metadata={"ffn_internal": True})
    dispersedDiameterModel: DispersedDiameterModel | None = None
    thermoResidualAlpha: int | float | None = None
    writeRestartFields: bool | None = None
