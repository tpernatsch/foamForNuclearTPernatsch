import numpy as np
from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import OffbeatDict, List
from foamForNuclear._attrs_tools import offbeat_define, _to_List_str
from typing import ClassVar
from attrs import field, validators as v


@offbeat_define
class DragModel(OffbeatDict):
    """
    Base class for fluid drag.
    """
    TYPE: ClassVar[str] = "none"
    zones: list | List = field(factory=list, metadata={"ffn_internal": True}, converter=_to_List_str)

    def __repr__(self, depth = 0):
        textZones = ""
        if (len(self.zones) > 0):
            textZones = "\"" + ':'.join(self.zones) + "\""

        return textZones + super().__repr__(depth)

class DragByRegime(DragModel):
    def __init__(
            self,
            regimeMap: str,
            zones=[]
        ):
        super().__init__("byRegime", zones)

        self.regimeMap = regimeMap
        self.regimes = []


    def __repr__(self, depth=0):
        self.__setitem__("regimeMap", self.regimeMap)

        for regimeName, dragModel in self.regimes:
            self.__setitem__(regimeName, dragModel)

        return super().__repr__(depth)

    def add_regime(self, dragModel: DragModel):
        check_type("dragModel", dragModel, DragModel)
        regimeName = dragModel.zones[0]
        dragModel.zones = []
        self.regimes.append((regimeName, dragModel))


@offbeat_define
class BaxiDalleDonne(DragModel):
    """
    Baxi Dalle Donne Drag Model

    Parameters
    ----------
    pinDiameter : int | float
    wireDiameter : int | float
    wireLeadLen : int | float
    """
    TYPE: ClassVar[str] = "BaxiDalleDonne"
    pinDiameter: int | float
    wireDiameter: int | float
    wireLeadLen: int | float


@offbeat_define
class Churchill(DragModel):
    """
    Churchill Drag Model

    Parameters
    ----------
    surfaceRoughness : int | float
    """
    TYPE: ClassVar[str] = "Churchill"
    surfaceRoughness: int | float


@offbeat_define
class Colebrook(DragModel):
    """
    Colebrook Drag Model:

    .. math::
        (\text{coeff} \times \log10(\text{Re}) + \text{const})^\text{exp}

    Parameters
    ----------
    coeff : int | float
    exp : int | float
    const : int | float
    """
    TYPE: ClassVar[str] = "Colebrook"
    coeff: int | float
    exp: int | float
    const: int | float

    def value(self, Re: float):
        return(pow(self.coeff * np.log10(Re) + self.const, self.exp))


@offbeat_define
class Engel(DragModel):
    """
    Engel Drag Model
    """
    TYPE: ClassVar[str] = "Engel"
    coeff: int | float


@offbeat_define
class ModifiedEngel(DragModel):
    """
    Modified Engel Drag Model
    """
    TYPE: ClassVar[str] = "ModifiedEngel"


@offbeat_define
class NoKazimiFluidStructure(DragModel):
    """
    No Kazimi Drag Model

    Parameters
    ----------
    pinDiameter : int | float
    wireDiameter : int | float
    wireLeadLen : int | float
    """
    TYPE: ClassVar[str] = "NoKazimi"
    pinDiameter: int | float
    wireDiameter: int | float
    wireLeadLen: int | float


@offbeat_define
class Rehme(DragModel):
    """
    Rehme Drag Model

    Parameters
    ----------
    numberOfPins : int
    pinDiameter : int | float
    wireDiameter : int | float
    wireLeadLen : int | float
    wetWrapPerimeter : int | float
    """
    TYPE: ClassVar[str] = "Rehme"
    numberOfPins: int | float
    pinDiameter: int | float
    wireDiameter: int | float
    wireLeadLen: int | float
    wetWrapPerimeter: int | float


@offbeat_define
class ReynoldsPower(DragModel):
    """
    Reynolds Power drag model.

    .. math::
        \text{coeff} \times \text{Re}^\text{exp} + \text{const}

    Parameters
    ----------
    coeff : int | float
    exp : int | float
    const : int | float, optional
        Default to 0
    """
    TYPE: ClassVar[str] = "ReynoldsPower"
    coeff: int | float
    exp: int | float
    const: int | float = 0.0

    def value(self, Re: float):
        return(self.coeff * pow(Re, self.exp) + self.const)


@offbeat_define
class Autruffe(DragModel):
    """
    Autruffe correlation for the friction factor between two fluids,
    originally developed for sodium liquid-vapour pairs. The Autruffe model only
    works for liquid-gas systems.

    See Autruffe, M.I., Wilson, G.J., Stewart, B., Kazimi, M.S., 1979. A proposed
    momentum exchange coefficient for two-phase modeling of sodium boiling.
    In: Proc. Int. Meeting on Fast Reactor Safety Tech. 4, Seattle, WA, pp.
    2515-2521.
    """
    TYPE: ClassVar[str] = "Autruffe"


@offbeat_define
class Bestion(DragModel):
    """
    Bestion Model for interfacial friction in fluid-fluid interactions.
    The Bestion model only works for liquid-gas systems.
    """
    TYPE: ClassVar[str] = "Bestion"


@offbeat_define
class BestionTRACE(DragModel):
    """
    Bestion model for interfacial friction in fluid-fluid interactions.
    Distribution parameter set to zero as in TRACE.
    The BestionTRACE model only works for liquid-gas systems.
    """
    TYPE: ClassVar[str] = "BestionTRACE"


@offbeat_define
class NoKazimiFluidFluid(DragModel):
    """
    No-Kazimi model for interfacial friction in fluid-fluid interactions.
    The NoKazimi model only works for liquid-gas systems.
    """
    TYPE: ClassVar[str] = "NoKazimi"
    pinPitch: int | float
    pinDiameter: int | float


@offbeat_define
class SchillerNaumann(DragModel):
    """
    Schiller-Neumann model for interfacial friction in fluid-fluid interactions.
    """
    TYPE: ClassVar[str] = "SchillerNaumann"

    def value(self, Re: float) -> float:
        if (Re < 1000):
            return(24 * (1.0 + 0.15*pow(Re, 0.687)) / Re)

        return(0.44)


@offbeat_define
class Wallis(DragModel):
    """
    Wallis correlation for interfacial liquid-vapour sodium friction in
    fluid-fluid interactions.
    The Wallis model only works for liquid-gas systems.

    H. Ninokata, T. Okano, "SABENA: subassembly boiling evolution numerical
    analysis", Nucl. Eng. Des., 120 (1990), pp. 349-367
    """
    TYPE: ClassVar[str] = "Wallis"