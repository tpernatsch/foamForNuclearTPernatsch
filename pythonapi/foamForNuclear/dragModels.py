import numpy as np
from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import OpenFOAMDict


_DRAG_MODEL_TYPES = {
    "BaxiDalleDonne", "Churchill", "Colebrook", "Engel", "NoKazimi", "Rehme",
    "ReynoldsPower", "byRegime", "constant", "modifiedEngel",
    "byRegime",
    "Autruffe", "Bestion", "BestionTRACE", "SchillerNaumann", "Wallis"
}


class DragModel(OpenFOAMDict):
    """
    Base class for fluid drag.
    """
    def __init__(
            self,
            type,
            zones=[],
        ):
        super().__init__()
        self.type = type
        self.zones = zones


    def __repr__(self, depth = 0):
        textZones = ""
        if (len(self.zones) > 0):
            textZones = "\"" + ':'.join(self.zones) + "\""

        return textZones + super().__repr__(depth)


    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _DRAG_MODEL_TYPES)
        self._type_ = type
        self.__setitem__("type", type)


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


class BaxiDalleDonne(DragModel):
    """
    Baxi Dalle Donne Drag Model

    Parameters
    ----------
    pinDiameter : float
    wireDiameter : float
    wireLeadLen : float
    """

    def __init__(
            self,
            pinDiameter: float,
            wireDiameter: float,
            wireLeadLen: float,
            zones: list[str]=[],
        ):
        super().__init__("BaxiDalleDonne", zones)

        self.pinDiameter = pinDiameter
        self.wireDiameter = wireDiameter
        self.wireLeadLen = wireLeadLen


    def __repr__(self, depth=0):
        self.__setitem__("pinDiameter", self.pinDiameter)
        self.__setitem__("wireDiameter", self.wireDiameter)
        self.__setitem__("wireLeadLen", self.wireLeadLen)

        return super().__repr__(depth)


class Churchill(DragModel):
    """
    Churchill Drag Model

    Parameters
    ----------
    surfaceRoughness : float
    """

    def __init__(
            self,
            surfaceRoughness: float,
            zones=[],
        ):
        super().__init__("Churchill", zones)

        self.surfaceRoughness = surfaceRoughness


    def __repr__(self, depth=0):
        self.__setitem__("surfaceRoughness", self.surfaceRoughness)

        return super().__repr__(depth)


class Colebrook(DragModel):
    r"""
    Colebrook Drag Model:

    .. math::
        (\text{coeff} \times \log10(\text{Re}) + \text{const})^\text{exp}

    Parameters
    ----------
    coeff : float
    exp : float
    const : float
    """

    def __init__(
            self,
            coeff: float,
            exp: float,
            const: float,
            zones=[],
        ):
        super().__init__("Colebrook", zones)

        self.coeff = coeff
        self.exp = exp
        self.const = const


    def __repr__(self, depth=0):
        self.__setitem__("coeff", self.coeff)
        self.__setitem__("exp", self.exp)
        self.__setitem__("const", self.const)

        return super().__repr__(depth)


    def value(self, Re: float):
        return(pow(self.coeff * np.log10(Re) + self.const, self.exp))


class Engel(DragModel):
    """
    Engel Drag Model
    """

    def __init__(
            self,
            zones=[],
        ):
        super().__init__("Engel", zones)


class ModifiedEngel(DragModel):
    """
    Modified Engel Drag Model
    """

    def __init__(
            self,
            zones=[],
        ):
        super().__init__("ModifiedEngel", zones)


class NoKazimiFluidStructureDragModel(DragModel):
    """
    No Kazimi Drag Model

    Parameters
    ----------
    pinDiameter : float
    wireDiameter : float
    wireLeadLen : float
    """

    def __init__(
            self,
            pinDiameter: float,
            wireDiameter: float,
            wireLeadLen: float,
            zones=[],
        ):
        super().__init__("NoKazimi", zones)

        self.pinDiameter = pinDiameter
        self.wireDiameter = wireDiameter
        self.wireLeadLen = wireLeadLen


    def __repr__(self, depth=0):
        self.__setitem__("pinDiameter", self.pinDiameter)
        self.__setitem__("wireDiameter", self.wireDiameter)
        self.__setitem__("wireLeadLen", self.wireLeadLen)

        return super().__repr__(depth)


class Rehme(DragModel):
    """
    Rehme Drag Model

    Parameters
    ----------
    numberOfPins : int
    pinDiameter : float
    wireDiameter : float
    wireLeadLen : float
    wetWrapPerimeter : float
    """

    def __init__(
            self,
            numberOfPins: int,
            pinDiameter: float,
            wireDiameter: float,
            wireLeadLen: float,
            wetWrapPerimeter: float,
            zones=[],
        ):
        super().__init__("Rehme", zones)

        self.numberOfPins = numberOfPins
        self.pinDiameter = pinDiameter
        self.wireDiameter = wireDiameter
        self.wireLeadLen = wireLeadLen
        self.wetWrapPerimeter = wetWrapPerimeter


    def __repr__(self, depth=0):
        self.__setitem__("numberOfPins", self.numberOfPins)
        self.__setitem__("pinDiameter", self.pinDiameter)
        self.__setitem__("wireDiameter", self.wireDiameter)
        self.__setitem__("wireLeadLen", self.wireLeadLen)
        self.__setitem__("wetWrapPerimeter", self.wetWrapPerimeter)

        return super().__repr__(depth)


class ReynoldsPower(DragModel):
    r"""
    Reynolds Power drag model.

    .. math::
        \text{coeff} \times \text{Re}^\text{exp} + \text{const}

    Parameters
    ----------
    coeff : float
    exp : float
    const : float, optional
        Default to 0
    """

    def __init__(
            self,
            coeff: float,
            exp: float,
            const: float=0,
            zones=[],
        ):
        super().__init__("ReynoldsPower", zones)

        self.coeff = coeff
        self.exp = exp
        self.const = const


    def __repr__(self, depth=0):
        self.__setitem__("coeff", self.coeff)
        self.__setitem__("exp", self.exp)
        self.__setitem__("const", self.const)

        return super().__repr__(depth)


    def value(self, Re: float):
        return(self.coeff * pow(Re, self.exp) + self.const)


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
    def __init__(self, zones=[]):
        super().__init__("Autruffe", zones)


class Bestion(DragModel):
    """
    Bestion Model for interfacial friction in fluid-fluid interactions.
    The Bestion model only works for liquid-gas systems.
    """
    def __init__(self, zones=[]):
        super().__init__("Bestion", zones)


class BestionTRACE(DragModel):
    """
    Bestion model for interfacial friction in fluid-fluid interactions.
    Distribution parameter set to zero as in TRACE.
    The BestionTRACE model only works for liquid-gas systems.
    """
    def __init__(self, zones=[]):
        super().__init__("BestionTRACE", zones)


class NoKazimiFluidFluidDragModel(DragModel):
    """
    No-Kazimi model for interfacial friction in fluid-fluid interactions.
    The NoKazimi model only works for liquid-gas systems.
    """
    def __init__(
            self,
            pinPitch: float,
            pinDiameter: float,
            zones=[]
        ):
        super().__init__("NoKazimi", zones)

        self.pinPitch = pinPitch
        self.pinDiameter = pinDiameter


    def __repr__(self, depth=0):
        self.__setitem__("pinPitch", self.pinPitch)
        self.__setitem__("pinDiameter", self.pinDiameter)

        return super().__repr__(depth)


class SchillerNaumann(DragModel):
    """
    Schiller-Neumann model for interfacial friction in fluid-fluid interactions.
    """
    def __init__(self, zones=[]):
        super().__init__("SchillerNaumann", zones)


    def value(self, Re: float) -> float:
        if (Re < 1000):
            return(24 * (1.0 + 0.15*pow(Re, 0.687)) / Re)

        return(0.44)


class Wallis(DragModel):
    """
    Wallis correlation for interfacial liquid-vapour sodium friction in
    fluid-fluid interactions.
    The Wallis model only works for liquid-gas systems.

    H. Ninokata, T. Okano, "SABENA: subassembly boiling evolution numerical
    analysis", Nucl. Eng. Des., 120 (1990), pp. 349-367
    """
    def __init__(self, zones=[]):
        super().__init__("Wallis", zones)
