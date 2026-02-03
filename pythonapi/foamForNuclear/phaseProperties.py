import numpy as np
from foamForNuclear.checkvalue import check_positive, check_type, check_value
from foamForNuclear.common import *
from foamForNuclear.dispersedDiameterModels import DispersedDiameterModel
from foamForNuclear.dragModels import DragModel
from foamForNuclear.heatExchangerModels import HeatExchangerModel
from foamForNuclear.heatTransferModels import HeatTransferModel
from foamForNuclear.openfoamFile import OpenFOAMFile
from foamForNuclear.pairGeometryModels import PairGeometryModel
from foamForNuclear.phaseChangeModel import PhaseChangeModel
from foamForNuclear.powerModels import PowerModel
from foamForNuclear.powerOffCriterionModels import PowerOffCriterionModel
from foamForNuclear.regimeMapModels import RegimeMapModel
from foamForNuclear.twoPhaseDragMultiplierModels import TwoPhaseDragMultiplierModel
from foamForNuclear.pump import Pump

_LATTICE_TYPES = {"square", "hexagon"}
_STATE_OF_MATTER_TYPES = {"liquid", "gas"}



class PassiveProperties(OpenFOAMDict):
    def __init__(
            self,
            volumetricArea,
            T,
            rho=None,
            Cp=None,
            rhoCp=None,
        ):
        super().__init__()
        self.volumetricArea = volumetricArea
        self.rho = rho
        self.Cp = Cp
        self.rhoCp = rhoCp
        self.T = T

    def __repr__(self, depth = 0):
        self.__setitem__('volumetricArea', self.volumetricArea)
        if (self.rho is not None and self.Cp is not None and self.rhoCp is not None):
            msg = "In PassiveProperties, provide only rho and Cp, OR rhoCp"
            raise ValueError(msg)
        if (self.rho is not None and self.Cp is not None):
            self.__setitem__('rho', self.rho)
            self.__setitem__('Cp', self.Cp)
        elif (self.rhoCp is not None):
            self.__setitem__('rhoCp', self.rhoCp)
        else:
            msg = "In PassiveProperties, provide either rho and Cp, OR rhoCp"
            raise ValueError(msg)
        self.__setitem__('T', self.T)

        return super().__repr__(depth)


class StructureProperty(OpenFOAMDict):
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
    def __init__(
            self,
            zones=[],
            volumeFraction: float=None,
            Dh: float=None,
            pitch: float=None,
            elementDiameter: float=None,
            latticeType: str=None,
            spacerThickness: float=0,
            gridFraction: float=0,
            wireDiameter: float=0,
            localX: Vector=None,
            localZ: Vector=None,
            localDhAnisotropy: Vector=None,
            localTortuosity: Vector=None,
        ):
        super().__init__()
        self.zones = zones
        self.volumeFraction = volumeFraction
        self.Dh = Dh
        self.localX = localX
        self.localZ = localZ
        self.localDhAnisotropy = localDhAnisotropy
        self.localTortuosity = localTortuosity

        self.pitch = pitch
        self.elementDiameter = elementDiameter
        self.latticeType = latticeType
        self.spacerThickness = spacerThickness
        self.gridFraction = gridFraction
        self.wireDiameter = wireDiameter

        self.powerModel = None
        self.passiveProperties = None
        self.pump = None


    def __repr__(self, depth = 0):
        textZones = "\"" + ':'.join(self.zones) + "\""

        if (self.powerModel is not None):
            self.__setitem__('powerModel', self.powerModel)

        if (self.passiveProperties is not None):
            self.__setitem__('passiveProperties', self.passiveProperties)

        if (self.pump is not None):
            self.__setitem__('pump', self.pump)

        return textZones + super().__repr__(depth)


    @property
    def zones(self):
        return self._zones

    @zones.setter
    def zones(self, zones) -> None:
        if (zones is not None):
            check_type("zones", zones, list)
            self._zones = zones
        else:
            self._zones = []

    @property
    def volumeFraction(self):
        return self._volumeFraction

    @volumeFraction.setter
    def volumeFraction(self, volumeFraction) -> None:
        check_type("volumeFraction", volumeFraction, (float, int), none_ok=True)
        if (volumeFraction is not None and (volumeFraction < 0 or 1 < volumeFraction)):
            msg = "volumeFraction must be between 0 and 1"
            raise ValueError(msg)
        self._volumeFraction = volumeFraction
        self.__setitem__('volumeFraction', volumeFraction)

    @property
    def Dh(self):
        return self._Dh

    @Dh.setter
    def Dh(self, Dh) -> None:
        check_type("Dh", Dh, (float, int), none_ok=True)
        if (Dh is not None):
            check_positive("Dh", Dh)
        self._Dh = Dh
        self.__setitem__('Dh', Dh)

    @property
    def pitch(self):
        return self._pitch

    @pitch.setter
    def pitch(self, pitch) -> None:
        check_type("pitch", pitch, (float, int), none_ok=True)
        if (pitch is not None):
            check_positive("pitch", pitch)
        self._pitch = pitch
        try:
            self.compute_hydraulic_parameters(
                self.pitch, self.elementDiameter, self.latticeType,
                self.spacerThickness, self.gridFraction, self.wireDiameter
            )
        except:
            pass

    @property
    def elementDiameter(self):
        return self._elementDiameter

    @elementDiameter.setter
    def elementDiameter(self, elementDiameter) -> None:
        check_type("elementDiameter", elementDiameter, (float, int), none_ok=True)
        if (elementDiameter is not None):
            check_positive("elementDiameter", elementDiameter)
        self._elementDiameter = elementDiameter
        try:
            self.compute_hydraulic_parameters(
                self.pitch, self.elementDiameter, self.latticeType,
                self.spacerThickness, self.gridFraction, self.wireDiameter
            )
        except:
            pass

    @property
    def latticeType(self):
        return self._latticeType

    @latticeType.setter
    def latticeType(self, latticeType) -> None:
        check_type("latticeType", latticeType, str, none_ok=True)
        if (latticeType is not None):
            check_value("latticeType", latticeType, _LATTICE_TYPES)
        self._latticeType = latticeType
        try:
            self.compute_hydraulic_parameters(
                self.pitch, self.elementDiameter, self.latticeType,
                self.spacerThickness, self.gridFraction, self.wireDiameter
            )
        except:
            pass

    @property
    def spacerThickness(self):
        return self._spacerThickness

    @spacerThickness.setter
    def spacerThickness(self, spacerThickness) -> None:
        check_type("spacerThickness", spacerThickness, (float, int), none_ok=True)
        if (spacerThickness is not None):
            check_positive("spacerThickness", spacerThickness)
        self._spacerThickness = spacerThickness
        try:
            self.compute_hydraulic_parameters(
                self.pitch, self.elementDiameter, self.latticeType,
                self.spacerThickness, self.gridFraction, self.wireDiameter
            )
        except:
            pass

    @property
    def gridFraction(self):
        return self._gridFraction

    @gridFraction.setter
    def gridFraction(self, gridFraction) -> None:
        check_type("gridFraction", gridFraction, (float, int), none_ok=True)
        if (gridFraction is not None and (gridFraction < 0 or 1 < gridFraction)):
            msg = "gridFraction must be between 0 and 1"
            raise ValueError(msg)
        self._gridFraction = gridFraction
        try:
            self.compute_hydraulic_parameters(
                self.pitch, self.elementDiameter, self.latticeType,
                self.spacerThickness, self.gridFraction, self.wireDiameter
            )
        except:
            pass

    @property
    def wireDiameter(self):
        return self._wireDiameter

    @wireDiameter.setter
    def wireDiameter(self, wireDiameter) -> None:
        check_type("wireDiameter", wireDiameter, (float, int), none_ok=True)
        if (wireDiameter is not None):
            check_positive("wireDiameter", wireDiameter)
        self._wireDiameter = wireDiameter
        try:
            self.compute_hydraulic_parameters(
                self.pitch, self.elementDiameter, self.latticeType,
                self.spacerThickness, self.gridFraction, self.wireDiameter
            )
        except:
            pass

    @property
    def localX(self):
        return self._localX

    @localX.setter
    def localX(self, localX) -> None:
        check_type("localX", localX, Vector, none_ok=True)
        self._localX = localX
        if (localX is not None):
            self.__setitem__('localX', localX)

    @property
    def localZ(self):
        return self._localZ

    @localZ.setter
    def localZ(self, localZ) -> None:
        check_type("localZ", localZ, Vector, none_ok=True)
        self._localZ = localZ
        if (localZ is not None):
            self.__setitem__('localZ', localZ)

    @property
    def localDhAnisotropy(self):
        return self._localDhAnisotropy

    @localDhAnisotropy.setter
    def localDhAnisotropy(self, localDhAnisotropy) -> None:
        check_type("localDhAnisotropy", localDhAnisotropy, Vector, none_ok=True)
        self._localDhAnisotropy = localDhAnisotropy
        if (localDhAnisotropy is not None):
            self.__setitem__('localDhAnisotropy', localDhAnisotropy)

    @property
    def localTortuosity(self):
        return self._localTortuosity

    @localTortuosity.setter
    def localTortuosity(self, localTortuosity) -> None:
        check_type("localTortuosity", localTortuosity, Vector, none_ok=True)
        self._localTortuosity = localTortuosity
        if (localTortuosity is not None):
            self.__setitem__('localTortuosity', localTortuosity)

    @property
    def powerModel(self):
        return self._powerModel

    @powerModel.setter
    def powerModel(self, powerModel) -> None:
        check_type("powerModel", powerModel, PowerModel, none_ok=True)
        self._powerModel = powerModel

    @property
    def passiveProperties(self):
        return self._passiveProperties

    @passiveProperties.setter
    def passiveProperties(self, passiveProperties) -> None:
        check_type("passiveProperties", passiveProperties, PassiveProperties, none_ok=True)
        self._passiveProperties = passiveProperties

    @property
    def pump(self):
        return self._pump

    @pump.setter
    def pump(self, pump) -> None:
        check_type("pump", pump, Pump, none_ok=True)
        self._pump = pump


    def add_power_model(self, powerModel: PowerModel):
        self.powerModel = powerModel


    def add_passive_structure(self, volumetricArea, T, rho=None, Cp=None, rhoCp=None):
        self.passiveProperties = PassiveProperties(
            volumetricArea=volumetricArea,
            rho=rho,
            Cp=Cp,
            T=T,
            rhoCp=rhoCp
        )


    def add_pump(self, pump: Pump):
        self.pump = pump


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


class StructureProperties(OpenFOAMListDict):
    """
    Collect all structure model in the structureProperties sub-dict of the
    GeN-Foam phaseProperties file.

    Attributes
    ----------
    heatExchangerModels : OpenFOAMListDict
        List of heat exchanger models
    powerOffCriterionModel : PowerOffCriterionModel
        Power-off criterion
    """
    def __init__(self, items=None):
        super().__init__((StructureProperty, PowerOffCriterionModel, OpenFOAMListDict), "structureProperties", items)

        self.heatExchangerModels = OpenFOAMListDict(HeatExchangerModel, "heatExchangers")
        self.append(self.heatExchangerModels)

        self.powerOffCriterionModel = None

    @property
    def powerOffCriterionModel(self):
        return self._powerOffCriterionModel

    @powerOffCriterionModel.setter
    def powerOffCriterionModel(self, powerOffCriterionModel) -> None:
        check_type("powerOffCriterionModel", powerOffCriterionModel, PowerOffCriterionModel, none_ok=True)
        self._powerOffCriterionModel = powerOffCriterionModel
        if (powerOffCriterionModel is not None):
            self.append(powerOffCriterionModel)


    def append(
            self,
            item: StructureProperty | HeatExchangerModel | PowerOffCriterionModel | OpenFOAMListDict
        ):
        # Add once a OpenFOAMListDict for heatExchangerModels
        if (
            isinstance(item, OpenFOAMListDict)
            and any([type(e) == OpenFOAMListDict for e in self])
        ):
            msg = "Can only declare 1 OpenFOAMListDict which is restricted to heatExchangerModels"
            raise ValueError(msg)

        # Add once a power off criterion
        if (
            isinstance(item, PowerOffCriterionModel)
            and any([type(e) == PowerOffCriterionModel for e in self])
        ):
            return

        if (isinstance(item, (StructureProperty, PowerOffCriterionModel, OpenFOAMListDict))):
            super().append(item)
        elif (isinstance(item, HeatExchangerModel)):
            self.heatExchangerModels.append(item)


class FluidProperty(OpenFOAMDict):
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
    def __init__(
            self,
            stateOfMatter: str,
            thermoResidualAlpha: float=0,
            residualAlpha: float=1e-9,
            writeRestartFields: bool=False,
            dispersedDiameterModel: DispersedDiameterModel=None,
        ):
        super().__init__()

        self.stateOfMatter = stateOfMatter
        self.thermoResidualAlpha = thermoResidualAlpha
        self.residualAlpha = residualAlpha
        self.writeRestartFields = writeRestartFields
        self.dispersedDiameterModel = dispersedDiameterModel


    @property
    def thermoResidualAlpha(self):
        return self._thermoResidualAlpha

    @thermoResidualAlpha.setter
    def thermoResidualAlpha(self, thermoResidualAlpha) -> None:
        check_type("thermoResidualAlpha", thermoResidualAlpha, (float, int))
        self._thermoResidualAlpha = thermoResidualAlpha
        self.__setitem__('thermoResidualAlpha', thermoResidualAlpha)

    @property
    def stateOfMatter(self):
        return self._stateOfMatter

    @stateOfMatter.setter
    def stateOfMatter(self, stateOfMatter) -> None:
        check_type("stateOfMatter", stateOfMatter, str)
        check_value("stateOfMatter", stateOfMatter, _STATE_OF_MATTER_TYPES)
        self._stateOfMatter = stateOfMatter
        self.__setitem__('stateOfMatter', stateOfMatter)

    @property
    def residualAlpha(self):
        return self._residualAlpha

    @residualAlpha.setter
    def residualAlpha(self, residualAlpha) -> None:
        check_type("residualAlpha", residualAlpha, (float, int))
        self._residualAlpha = residualAlpha
        self.__setitem__('residualAlpha', residualAlpha)

    @property
    def writeRestartFields(self):
        return self._writeRestartFields

    @writeRestartFields.setter
    def writeRestartFields(self, writeRestartFields) -> None:
        check_type("writeRestartFields", writeRestartFields, bool)
        self._writeRestartFields = writeRestartFields
        self.__setitem__('writeRestartFields', writeRestartFields)

    @property
    def dispersedDiameterModel(self):
        return self._dispersedDiameterModel

    @dispersedDiameterModel.setter
    def dispersedDiameterModel(self, dispersedDiameterModel) -> None:
        check_type("dispersedDiameterModel", dispersedDiameterModel, DispersedDiameterModel, none_ok=True)
        self._dispersedDiameterModel = dispersedDiameterModel
        self.__setitem__('dispersedDiameterModel', dispersedDiameterModel)


class PhaseModels:
    def __init__(
            self,
            expected_type,
            name: str,
            phaseNames: list[str]=['fluid']
        ):
        self.expected_type = expected_type
        self.name = name

        self.modelsFluidStructure: dict[list] = {phaseName: [] for phaseName in phaseNames}

        self.phaseNames = phaseNames
        self.modelsFluidFluid: list = []


    def __repr__(self, depth = 0):
        # One phase solver
        if (len(self.phaseNames) == 1):
            obj = OpenFOAMListDict(
                self.expected_type,
                name=self.name,
                items=self.modelsFluidStructure[self.phaseNames[0]]
            )

            return(obj.__repr__(depth=depth))

        # Two phase solver
        obj = OpenFOAMListDict(
            (OpenFOAMListDict, self.expected_type),
            name=self.name
        )

        # Fluid-fluid
        if (len(self.modelsFluidFluid) > 0):
            # I don't like this distinction, the dict in phaseProperties should
            # be more uniform
            if (self.expected_type == DragModel):
                for model in self.modelsFluidFluid:
                    obj.append(model)
            else:
                obj.append(OpenFOAMListDict(
                    self.expected_type,
                    name='"'+'.'.join(self.phaseNames)+'"',
                    items=self.modelsFluidFluid
                ))

        # Fluid structure
        for phaseName, models in self.modelsFluidStructure.items():
            if (len(models) > 0):
                obj.append(OpenFOAMListDict(
                    self.expected_type,
                    name=f"\"{phaseName}.structure\"",
                    items=models
                ))

        return obj.__repr__(depth)

    @property
    def is_empty(self) -> bool:
        return(
            len(self.modelsFluidFluid) == 0
            and len(list(self.modelsFluidStructure.keys())) == 0
        )

    @property
    def phaseNames(self):
        return self._phaseNames

    @phaseNames.setter
    def phaseNames(self, phaseNames) -> None:
        check_type("phaseNames", phaseNames, list, none_ok=True)
        self._phaseNames = phaseNames
        if (self.phaseNames is None):
            self._phaseNames = ['fluid']

        self.modelsFluidStructure = {phaseName: [] for phaseName in self._phaseNames}


    def append(
            self,
            model,
            interactionType: str="fluid-structure",
            phaseName: str=None
        ):
        """
        Parameters
        ----------
        model
            Model to append
        interactionType : str
            Options: `"fluid-fluid"`, `"fluid-structure"`
            (default `fluid-structure`).
        phaseName : str
            Name of the phase used on the `twoPhase` solver (default `None`).
            Keep to `None` if fluid-fluid model.
        """
        check_type("model", model, self.expected_type)
        check_type("interactionType", interactionType, str)
        check_value("interactionType", interactionType, {"fluid-fluid", "fluid-structure"})
        check_type("phaseName", phaseName, str, none_ok=True)

        if (interactionType == "fluid-structure" and phaseName is not None):
            check_value("phaseName", phaseName, self.phaseNames)
            self.modelsFluidStructure[phaseName].append(model)

        elif (interactionType == "fluid-fluid"):
            self.modelsFluidFluid.append(model)

        else:
            self.modelsFluidStructure["fluid"].append(model)



class HeatTransferModels(PhaseModels):
    def __init__(
            self,
            phaseNames: list[str]=['fluid']
        ):
        super().__init__(
            expected_type=HeatTransferModel,
            name="heatTransferModels",
            phaseNames=phaseNames
        )


class DragModels(PhaseModels):
    def __init__(
            self,
            phaseNames: list[str]=['fluid']
        ):
        super().__init__(
            expected_type=DragModel,
            name="dragModels",
            phaseNames=phaseNames
        )


class PairGeometryModels(PhaseModels):
    def __init__(
            self,
            phaseNames: list[str]=['fluid']
        ):
        super().__init__(
            expected_type=PairGeometryModel,
            name="pairGeometryModels",
            phaseNames=phaseNames
        )


class PhaseProperties(OpenFOAMFile):
    """
    Phase properties object to assign the power, structure, drag, heat transfer,
    ... models.

    Parameters
    ----------
    structureProperties : StructureProperties
        (default `None`).
    regimeMapModels : OpenFOAMListDict
        (default `None`).
    dragModels : DragModels
        (default `None`).
    heatTransferModels : HeatTransferModels
        (default `None`).
    twoPhaseDragMultiplierModel : TwoPhaseDragMultiplierModel
        (default `None`).
    pairGeometryModels : PairGeometryModels
        (default `None`).
    phaseChangeModel: PhaseChangeModel
        (default `None`).
    residualKd : float
        Default `None`
    region : str
        Name of the region (default `""`).
    phaseNames : list[str]
        List of phase names used in two phase calculations (default `None`). If
        only one phase is at play in the calculations, this entry is not
        mandatory. Only used in `twoPhase` solver.
    fluidProperties : list[FluidProperty]
        Only used in `twoPhase` solver (default `[]`).
    virtualMassCoeff : float
        Only used in `twoPhase` solver (default `0.1`).
    """
    def __init__(
            self,
            structureProperties: StructureProperties=None,
            regimeMapModels=None,
            dragModels: DragModels=None,
            heatTransferModels: HeatTransferModels=None,
            twoPhaseDragMultiplierModel: TwoPhaseDragMultiplierModel=None,
            pairGeometryModels: PairGeometryModels=None,
            phaseChangeModel: PhaseChangeModel=None,
            residualKd: float=None,
            region: str="",
            phaseNames: list[str]=None,
            fluidProperties: list[FluidProperty]=[],
            virtualMassCoeff: float=0.1
        ):
        super().__init__("phaseProperties", "constant", region)

        # Common one and two phase
        self.structureProperties: StructureProperties = structureProperties
        self.regimeMapModels = regimeMapModels
        self.dragModels = dragModels
        self.heatTransferModels = heatTransferModels
        self.twoPhaseDragMultiplierModel = twoPhaseDragMultiplierModel
        self.pairGeometryModels: PairGeometryModels = pairGeometryModels
        self.phaseChangeModel = phaseChangeModel
        self.residualKd = residualKd

        # Two phase
        self.phaseNames = phaseNames
        self.fluidProperties = fluidProperties
        self.virtualMassCoeff = virtualMassCoeff


    def isTwoPhase(self) -> bool:
        return(self.phaseNames is not None)

    @property
    def residualKd(self):
        return self._residualKd

    @residualKd.setter
    def residualKd(self, residualKd) -> None:
        check_type("residualKd", residualKd, (float, int), none_ok=True)
        self._residualKd = residualKd

    @property
    def structureProperties(self):
        return self._structureProperties

    @structureProperties.setter
    def structureProperties(self, structureProperties) -> None:
        if structureProperties is not None:
            check_type("structureProperties", structureProperties, StructureProperties)
            self._structureProperties = structureProperties
        else:
            self._structureProperties = StructureProperties()

    @property
    def regimeMapModels(self):
        return self._regimeMapModels

    @regimeMapModels.setter
    def regimeMapModels(self, regimeMapModels) -> None:
        if regimeMapModels is not None:
            check_type("regimeMapModels", regimeMapModels, OpenFOAMListDict)
            self._regimeMapModels = regimeMapModels
        else:
            self._regimeMapModels = OpenFOAMListDict(RegimeMapModel, "regimeMapModels")

    @property
    def dragModels(self):
        return self._dragModels

    @dragModels.setter
    def dragModels(self, dragModels) -> None:
        if dragModels is not None:
            check_type("dragModels", dragModels, DragModels)
            self._dragModels = dragModels
        else:
            self._dragModels = DragModels()

    @property
    def heatTransferModels(self):
        return self._heatTransferModels

    @heatTransferModels.setter
    def heatTransferModels(self, heatTransferModels) -> None:
        if heatTransferModels is not None:
            check_type("heatTransferModels", heatTransferModels, HeatTransferModels)
            self._heatTransferModels = heatTransferModels
        else:
            self._heatTransferModels = HeatTransferModels()

    @property
    def twoPhaseDragMultiplierModel(self):
        return self._twoPhaseDragMultiplierModel

    @twoPhaseDragMultiplierModel.setter
    def twoPhaseDragMultiplierModel(self, twoPhaseDragMultiplierModel) -> None:
        check_type("twoPhaseDragMultiplierModel", twoPhaseDragMultiplierModel, TwoPhaseDragMultiplierModel, none_ok=True)
        self._twoPhaseDragMultiplierModel = twoPhaseDragMultiplierModel

    @property
    def pairGeometryModels(self):
        return self._pairGeometryModels

    @pairGeometryModels.setter
    def pairGeometryModels(self, pairGeometryModels) -> None:
        if pairGeometryModels is not None:
            check_type("pairGeometryModels", pairGeometryModels, PairGeometryModels)
            self._pairGeometryModels = pairGeometryModels
        else:
            self._pairGeometryModels = PairGeometryModels()

    @property
    def phaseChangeModel(self):
        return self._phaseChangeModel

    @phaseChangeModel.setter
    def phaseChangeModel(self, phaseChangeModel) -> None:
        check_type("phaseChangeModel", phaseChangeModel, PhaseChangeModel, none_ok=True)
        self._phaseChangeModel = phaseChangeModel

    @property
    def phaseNames(self):
        return self._phaseNames

    @phaseNames.setter
    def phaseNames(self, phaseNames) -> None:
        check_type("phaseNames", phaseNames, list, none_ok=True)
        self._phaseNames = phaseNames
        self.heatTransferModels.phaseNames = phaseNames
        self.dragModels.phaseNames = phaseNames
        self.pairGeometryModels.phaseNames = phaseNames

    @property
    def virtualMassCoeff(self):
        return self._virtualMassCoeff

    @virtualMassCoeff.setter
    def virtualMassCoeff(self, virtualMassCoeff) -> None:
        check_type("virtualMassCoeff", virtualMassCoeff, float)
        self._virtualMassCoeff = virtualMassCoeff


    def add_fluid_properties(self, fluidProperty: FluidProperty):
        check_type("fluidProperty", fluidProperty, FluidProperty)
        self.fluidProperties.append(fluidProperty)


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = ""

        if (self.isTwoPhase()):
            for i, phaseName in enumerate(self.phaseNames):
                text += addParameter(f"fluid{i+1}", phaseName, isAddExtraLine=True)

            if (len(self.fluidProperties) != len(self.phaseNames)):
                msg = "fluidProperties and phaseNames are different sizes. Add or remove fluid properties to use the two phase solver."
                raise ValueError(msg)

            for phaseName, fluidProperty in zip(self.phaseNames, self.fluidProperties):
                text += f"{phaseName}Properties{fluidProperty!r}\n"

            text += addParameter("virtualMassCoeff", self.virtualMassCoeff, isAddExtraLine=True)


        text += f"{self.structureProperties!r}\n"

        text += f"{self.regimeMapModels!r}\n"

        text += "physicsModels\n"
        text += "{\n"
        text += f"{self.dragModels.__repr__(depth=1)}\n"
        text += f"{self.heatTransferModels.__repr__(depth=1)}\n"
        if (self.twoPhaseDragMultiplierModel is not None):
            text += f"{tab}twoPhaseDragMultiplierModel{self.twoPhaseDragMultiplierModel.__repr__(depth=1)}\n"
        if (not self.pairGeometryModels.is_empty):
            text += f"{self.pairGeometryModels.__repr__(depth=1)}\n"
        if (self.phaseChangeModel is not None):
            text += f"{tab}{self.phaseChangeModel.__repr__(depth=1)}\n"
        text += "}\n\n"

        if (self.residualKd is not None):
            text += addParameter("residualKd", self.residualKd, isAddExtraLine=True)

        return(text)
