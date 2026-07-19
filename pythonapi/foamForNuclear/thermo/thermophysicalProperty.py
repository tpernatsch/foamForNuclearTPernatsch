from abc import abstractmethod
from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import OpenFOAMDict, addParameter
from foamForNuclear.openfoamFile import OpenFOAMFile


class BaseThermophysicalProperty(OpenFOAMFile):
    """
    Base thermophysical properties object.

    Parameters
    ----------
    region : str
        Name of the region
    ext : str
        Extension at the end of the file, e.g `".liquid"` (default `""`)

    Attributes
    ----------
    region : str
        Name of the region
    pRef : float
        Reference pressure
    """

    def __init__(self, region: str = "", ext: str = ""):
        super().__init__("thermophysicalProperties", "constant", region, ext=ext)

        if not hasattr(self, "thermoType"):
            self.thermoType = ThermoType()

        if not hasattr(self, "description"):
            self.description = None

        if not hasattr(self, "pRef"):
            self.pRef = 1e5
        if not hasattr(self, "molWeight"):
            self.molWeight = None
        if not hasattr(self, "rho"):
            self.rho = None
        if not hasattr(self, "Cp"):
            self.Cp = None
        if not hasattr(self, "Hf"):
            self.Hf = None
        if not hasattr(self, "Sf"):
            self.Sf = None
        if not hasattr(self, "mu"):
            self.mu = None
        if not hasattr(self, "kappa"):
            self.kappa = None
        if not hasattr(self, "Pr"):
            self.Pr = None

        if not hasattr(self, "mixture"):
            self.mixture = OpenFOAMDict()


    @abstractmethod
    def update_mixture(self):
        pass


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = ""
        text += f"thermoType{self.thermoType!r}\n"

        text += addParameter("pRef", self.pRef, isAddExtraLine=True)

        if (self.description is not None):
            for line in self.description.split("\n"):
                text += f"// {line}\n"

        self.update_mixture()
        text += f"mixture{self.mixture!r}\n"

        return(text)


_THERMO_TYPE_TYPE_TYPES = {"heRhoThermo"}
_THERMO_TYPE_MIXTURE_TYPES = {
    "homogeneousMixture", "inhomogeneousMixture", "multiComponentMixture",
    "pureMixture", "pureZoneMixture", "reactingMixture",
    "singleStepReactingMixture", "veryInhomogeneousMixture"
}
_THERMO_TYPE_TRANSPORT_TYPES = {"const", "sutherland", "polynomial", "WLF", "tabulated"}
_THERMO_TYPE_THERMO_TYPES = {"hConst", "janaf", "eConst", "hPolynomial", "hTabulated"}
_THERMO_TYPE_EQUATIONOFSTATE_TYPES = {
    "perfectGas", "incompressiblePerfectGas", "PengRobinsonGas",
    "adiabaticPerfectFluid", "perfectFluid", "rhoConst", "icoPolynomial",
    "rPolynomial", "Boussinesq", "icoTabulated"
}
_THERMO_TYPE_SPECIE_TYPES = {"specie"}
_THERMO_TYPE_ENERGY_TYPES = {"sensibleEnthalpy", "sensibleInternalEnergy"}
_THERMO_TYPE_PROPERTIES_TYPES = {"liquid"}


class ThermoType(OpenFOAMDict):
    """
    Thermophysical properties types.

    Parameters
    ----------
    type : str
        General thermophysical model calculation based on enthalpy h or internal
        energy e, and density rho (default `heRhoThermo`).
    mixture : str
        Mixture type (default `pureMixture`).
    transport : str
        Transport type.
    thermo : str
        Thermodynamics type.
    equationOfState : str
        Boussinesq `rho(T) = rho_0 (1 - beta(T-T_0))`
    specie : str
        Thermophysical properties of species, derived from Cp, h and/or s
        (default `specie`).
    energy : str
        Energy type (default `sensibleEnthalpy`).
    properties : str
        Optional Thermophysical properties implemented in the OpenFOAM source
        code format.
    """

    def __init__(
            self,
            type: str='heRhoThermo',
            mixture: str='pureMixture',
            transport: str=None,
            thermo: str=None,
            equationOfState: str=None,
            specie: str=None,
            energy: str='sensibleEnthalpy',
            properties: str=None,
        ):
        super().__init__()
        self.type = type
        self.mixture = mixture
        self.transport = transport
        self.thermo = thermo
        self.equationOfState = equationOfState
        self.specie = specie
        self.energy = energy
        self.properties = properties

        # super().__init__(name="thermoType")


    def __repr__(self, depth: int=0):
        if (self.mixture is not None):
            self.__setitem__('mixture', self.mixture)
        if (self.transport is not None):
            self.__setitem__('transport', self.transport)
        if (self.thermo is not None):
            self.__setitem__('thermo', self.thermo)
        if (self.equationOfState is not None):
            self.__setitem__('equationOfState', self.equationOfState)
        if (self.specie is not None):
            self.__setitem__('specie', self.specie)
        if (self.energy is not None):
            self.__setitem__('energy', self.energy)
        if (self.properties is not None):
            self.__setitem__('properties', self.properties)

        return super().__repr__(depth)


    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _THERMO_TYPE_TYPE_TYPES)
        self._type = type
        self.__setitem__('type', self.type)

    @property
    def mixture(self):
        return self._mixture

    @mixture.setter
    def mixture(self, mixture) -> None:
        check_type("mixture", mixture, str)
        check_value("mixture", mixture, _THERMO_TYPE_MIXTURE_TYPES)
        self._mixture = mixture

    @property
    def transport(self):
        return self._transport

    @transport.setter
    def transport(self, transport) -> None:
        check_type("transport", transport, str, none_ok=True)
        if transport is not None:
            check_value("transport", transport, _THERMO_TYPE_TRANSPORT_TYPES)
        self._transport = transport

    @property
    def thermo(self):
        return self._thermo

    @thermo.setter
    def thermo(self, thermo) -> None:
        check_type("thermo", thermo, str, none_ok=True)
        if thermo is not None:
            check_value("thermo", thermo, _THERMO_TYPE_THERMO_TYPES)
        self._thermo = thermo

    @property
    def equationOfState(self):
        return self._equationOfState

    @equationOfState.setter
    def equationOfState(self, equationOfState) -> None:
        check_type("equationOfState", equationOfState, str, none_ok=True)
        if equationOfState is not None:
            check_value("equationOfState", equationOfState, _THERMO_TYPE_EQUATIONOFSTATE_TYPES)
        self._equationOfState = equationOfState

    @property
    def specie(self):
        return self._specie

    @specie.setter
    def specie(self, specie) -> None:
        check_type("specie", specie, str, none_ok=True)
        if specie is not None:
            check_value("specie", specie, _THERMO_TYPE_SPECIE_TYPES)
        self._specie = specie

    @property
    def energy(self):
        return self._energy

    @energy.setter
    def energy(self, energy) -> None:
        check_type("energy", energy, str, none_ok=True)
        if energy is not None:
            check_value("energy", energy, _THERMO_TYPE_ENERGY_TYPES)
        self._energy = energy

    @property
    def properties(self):
        return self._properties

    @properties.setter
    def properties(self, properties) -> None:
        check_type("properties", properties, str, none_ok=True)
        if properties is not None:
            check_value("properties", properties, _THERMO_TYPE_PROPERTIES_TYPES)
        self._properties = properties
