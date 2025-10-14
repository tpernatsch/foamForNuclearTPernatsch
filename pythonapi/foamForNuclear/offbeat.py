import numpy as np
from foamForNuclear.axialProfileModels import AxialProfile
from foamForNuclear.azimuthalProfileModels import AzimuthalProfile
from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import *
from foamForNuclear.fvSchemes import fvSchemes
from foamForNuclear.fvSolution import fvSolution, fvSolutionSolver
from foamForNuclear.mesh.mesh import Mesh
from foamForNuclear.radialProfileModels import RadialProfile
from foamForNuclear.solver import Solver
from foamForNuclear.thermomechanicalMaterial.thermomechanicalMaterial import BaseThermomechanicalMaterial
from foamForNuclear.timeFolder import TimeFolder


_OFFBEAT_SOLVER_TYPES = {"offbeat", "extendedThermoMechanics"}

_OFFBEAT_THERMAL_SOLVER_TYPES = {
    "solidConduction", "readTemperature", "fromLatestTime"
}
_OFFBEAT_MECHANICAL_SOLVER_TYPES = {
    "largeStrainTotLag", "fromLatestTime", "largeStrainUpdLag", "smallStrain",
    "smallStrainIncrementalUpdated", "smallStrainIncremental"
}
_OFFBEAT_NEUTRONICS_SOLVER_TYPES = {"diffusion", "fromLatestTime"}
_OFFBEAT_ELEMENT_TRANSPORT_SOLVER_TYPES = {"byList", "fromLatestTime"}

_OFFBEAT_MATERIAL_PROPERTIES_TYPES = {"byZone"}
_OFFBEAT_RHEOLOGY_TYPES = {"byMaterial"}

_OFFBEAT_HEAT_SOURCE_TYPES = {
    "fmiLhgr", "constantLhgr", "timeDependentVhgr", "fromLatestTime",
    "timeDependentLhgr"
}
_OFFBEAT_BURNUP_TYPES = {"none", "fromPower", "fromLatestTime", "Lassmann"}
_OFFBEAT_FAST_FLUX_TYPES = {
    "none", "fromLatestTime", "timeDependentAxialProfile"
}
_OFFBEAT_CORROSION_TYPES = {"fromLatestTime"}
_OFFBEAT_GAP_GAS_TYPES = {"FRAPCON", "none", "TRISO"}
_OFFBEAT_FGR_TYPES = {"none", "SCIANTIX"}
_OFFBEAT_SLICE_MAPPER_TYPES = {
    "autoAxialSlices", "none", "byPellets", "byMaterial"
}

_INTERPOLATION_METHOD_TYPES = {"linear"}



class ThermoMechanicsCouplingOptions(OpenFOAMDict):
    def __init__(
            self,
            correctTFromTH: bool=None,
            correctDispForNeutro: bool=None
        ):
        super().__init__()
        self.correctTFromTH = correctTFromTH
        self.correctDispForNeutro = correctDispForNeutro


    @property
    def correctTFromTH(self):
        return self._correctTFromTH

    @correctTFromTH.setter
    def correctTFromTH(self, correctTFromTH) -> None:
        if (correctTFromTH is not None):
            check_type("correctTFromTH", correctTFromTH, bool)
            self._correctTFromTH = correctTFromTH
        else:
            self._correctTFromTH = True
        self.__setitem__('correctTFromTH', self.correctTFromTH)

    @property
    def correctDispForNeutro(self):
        return self._correctDispForNeutro

    @correctDispForNeutro.setter
    def correctDispForNeutro(self, correctDispForNeutro) -> None:
        if (correctDispForNeutro is not None):
            check_type("correctDispForNeutro", correctDispForNeutro, bool)
            self._correctDispForNeutro = correctDispForNeutro
        else:
            self._correctDispForNeutro = True
        self.__setitem__('correctDispForNeutro', self.correctDispForNeutro)


class GlobalOptions(OpenFOAMDict):
    def __init__(
            self,
            pinDirection: Vector=Vector(0, 0, 1),
            reactorType: str="LWR",
            angularFraction: float=0,
            linkedFuel: bool=False
        ):
        super().__init__()

        self.pinDirection = pinDirection
        self.reactorType = reactorType
        self.angularFraction = angularFraction
        self.linkedFuel = linkedFuel


    @property
    def pinDirection(self):
        return self._pinDirection

    @pinDirection.setter
    def pinDirection(self, pinDirection) -> None:
        check_type("pinDirection", pinDirection, Vector)
        self._pinDirection = pinDirection
        self.__setitem__('pinDirection', self.pinDirection)


    @property
    def reactorType(self):
        return self._reactorType

    @reactorType.setter
    def reactorType(self, reactorType) -> None:
        check_type("reactorType", reactorType, str)
        self._reactorType = reactorType
        self.__setitem__('reactorType', self.reactorType)


    @property
    def angularFraction(self):
        return self._angularFraction

    @angularFraction.setter
    def angularFraction(self, angularFraction) -> None:
        check_type("angularFraction", angularFraction, (float, int))
        self._angularFraction = angularFraction
        self.__setitem__('angularFraction', self.angularFraction)


    @property
    def linkedFuel(self):
        return self._linkedFuel

    @linkedFuel.setter
    def linkedFuel(self, linkedFuel) -> None:
        check_type("linkedFuel", linkedFuel, bool)
        self._linkedFuel = linkedFuel
        self.__setitem__('linkedFuel', self.linkedFuel)


#==============================================================================*
# Thermal solver

class ThermalSolverOptions(OpenFOAMDict):
    def __init__(
            self,
            type: str="fromLatestTime",
            heatFluxSummary: bool=False,
            calculateEnthalpy: bool=False,
            averageFieldsForMappings: bool=True,
            # patchForAverage: str=None,
        ):
        super().__init__()

        self.type = type
        self.heatFluxSummary = heatFluxSummary
        self.calculateEnthalpy = calculateEnthalpy
        self.averageFieldsForMappings = averageFieldsForMappings
        # self.patchForAverage = patchForAverage


    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _OFFBEAT_THERMAL_SOLVER_TYPES)
        self._type = type


    @property
    def heatFluxSummary(self):
        return self._heatFluxSummary

    @heatFluxSummary.setter
    def heatFluxSummary(self, heatFluxSummary) -> None:
        check_type("heatFluxSummary", heatFluxSummary, bool)
        self._heatFluxSummary = heatFluxSummary
        self.__setitem__('heatFluxSummary', self.heatFluxSummary)


    @property
    def calculateEnthalpy(self):
        return self._calculateEnthalpy

    @calculateEnthalpy.setter
    def calculateEnthalpy(self, calculateEnthalpy) -> None:
        check_type("calculateEnthalpy", calculateEnthalpy, bool)
        self._calculateEnthalpy = calculateEnthalpy
        self.__setitem__('calculateEnthalpy', self.calculateEnthalpy)


    @property
    def averageFieldsForMappings(self):
        return self._averageFieldsForMappings

    @averageFieldsForMappings.setter
    def averageFieldsForMappings(self, averageFieldsForMappings) -> None:
        check_type("averageFieldsForMappings", averageFieldsForMappings, bool)
        self._averageFieldsForMappings = averageFieldsForMappings
        self.__setitem__('averageFieldsForMappings', self.averageFieldsForMappings)


    # @property
    # def patchForAverage(self):
    #     return self._patchForAverage

    # @patchForAverage.setter
    # def patchForAverage(self, patchForAverage) -> None:
    #     check_type("patchForAverage", patchForAverage, str, none_ok=True)
    #     self._patchForAverage = patchForAverage
    #     if (patchForAverage is not None):
    #         self.__setitem__('patchForAverage', self.patchForAverage)


class ReadTemperatureThermalSolverOptions(ThermalSolverOptions):
    def __init__(
            self,
            heatFluxSummary: bool=False,
            calculateEnthalpy: bool=False,
            averageFieldsForMappings: bool=True):
        super().__init__(
            type="readTemperature",
            heatFluxSummary=heatFluxSummary,
            calculateEnthalpy=calculateEnthalpy,
            averageFieldsForMappings=averageFieldsForMappings
        )


class SolidConductionThermalSolverOptions(ThermalSolverOptions):
    def __init__(
            self,
            heatFluxSummary: bool=False,
            calculateEnthalpy: bool=False,
            averageFieldsForMappings: bool=True):
        super().__init__(
            type="solidConduction",
            heatFluxSummary=heatFluxSummary,
            calculateEnthalpy=calculateEnthalpy,
            averageFieldsForMappings=averageFieldsForMappings
        )


#==============================================================================*
# Mechanical solver

class MechanicsSolverOptions(OpenFOAMDict):
    def __init__(
            self,
            type: str="fromLatestTime",
            forceSummary: bool=False,
            cylindricalStress: bool=False,
            sphericalStress: bool=False,
            RhieChowCorrection: bool=True,
        ):
        super().__init__()

        self.type = type
        self.forceSummary = forceSummary
        self.cylindricalStress = cylindricalStress
        self.sphericalStress = sphericalStress
        self.RhieChowCorrection = RhieChowCorrection
        self.multiMaterialCorrection = OpenFOAMDict({
            "type": "uniform",
            "defaultWeights": 1
        })

    def __repr__(self, depth = 0):
        self.__setitem__("multiMaterialCorrection", self.multiMaterialCorrection)

        return super().__repr__(depth)


    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _OFFBEAT_MECHANICAL_SOLVER_TYPES)
        self._type = type


    @property
    def forceSummary(self):
        return self._forceSummary

    @forceSummary.setter
    def forceSummary(self, forceSummary) -> None:
        check_type("forceSummary", forceSummary, bool)
        self._forceSummary = forceSummary
        self.__setitem__('forceSummary', self.forceSummary)


    @property
    def cylindricalStress(self):
        return self._cylindricalStress

    @cylindricalStress.setter
    def cylindricalStress(self, cylindricalStress) -> None:
        check_type("cylindricalStress", cylindricalStress, bool)
        self._cylindricalStress = cylindricalStress
        self.__setitem__('cylindricalStress', self.cylindricalStress)


    @property
    def sphericalStress(self):
        return self._sphericalStress

    @sphericalStress.setter
    def sphericalStress(self, sphericalStress) -> None:
        check_type("sphericalStress", sphericalStress, bool)
        self._sphericalStress = sphericalStress
        self.__setitem__('sphericalStress', self.sphericalStress)


    @property
    def RhieChowCorrection(self):
        return self._RhieChowCorrection

    @RhieChowCorrection.setter
    def RhieChowCorrection(self, RhieChowCorrection) -> None:
        check_type("RhieChowCorrection", RhieChowCorrection, bool)
        self._RhieChowCorrection = RhieChowCorrection
        self.__setitem__('RhieChowCorrection', self.RhieChowCorrection)


class SmallStrainMechanicsSolverOptions(MechanicsSolverOptions):
    def __init__(
            self,
            forceSummary: bool=False,
            cylindricalStress: bool=False,
            sphericalStress: bool=False,
            RhieChowCorrection: bool=True,
            inertialRelaxation: bool=False,
            inertialRelaxationDeltaT: float=1e-7,
        ):
        super().__init__(
            type="smallStrain",
            forceSummary=forceSummary,
            cylindricalStress=cylindricalStress,
            sphericalStress=sphericalStress,
            RhieChowCorrection=RhieChowCorrection
        )

        self.inertialRelaxation = inertialRelaxation
        self.inertialRelaxationDeltaT = inertialRelaxationDeltaT


    @property
    def inertialRelaxation(self):
        return self._inertialRelaxation

    @inertialRelaxation.setter
    def inertialRelaxation(self, inertialRelaxation) -> None:
        check_type("inertialRelaxation", inertialRelaxation, bool)
        self._inertialRelaxation = inertialRelaxation
        self.__setitem__('inertialRelaxation', self.inertialRelaxation)


    @property
    def inertialRelaxationDeltaT(self):
        return self._inertialRelaxationDeltaT

    @inertialRelaxationDeltaT.setter
    def inertialRelaxationDeltaT(self, inertialRelaxationDeltaT) -> None:
        check_type("inertialRelaxationDeltaT", inertialRelaxationDeltaT, (float, int))
        self._inertialRelaxationDeltaT = inertialRelaxationDeltaT
        self.__setitem__('inertialRelaxationDeltaT', self.inertialRelaxationDeltaT)


class SmallStrainIncrementalUpdatedMechanicsSolverOptions(MechanicsSolverOptions):
    def __init__(
            self,
            forceSummary: bool=False,
            cylindricalStress: bool=False,
            sphericalStress: bool=False,
            RhieChowCorrection: bool=True,
            updateMesh: bool=True
        ):
        super().__init__(
            type="smallStrainIncrementalUpdated",
            forceSummary=forceSummary,
            cylindricalStress=cylindricalStress,
            sphericalStress=sphericalStress,
            RhieChowCorrection=RhieChowCorrection
        )

        self.updateMesh = updateMesh


    @property
    def updateMesh(self):
        return self._updateMesh

    @updateMesh.setter
    def updateMesh(self, updateMesh) -> None:
        check_type("updateMesh", updateMesh, bool)
        self._updateMesh = updateMesh
        self.__setitem__('updateMesh', self.updateMesh)


class LargeStrainTotLagMechanicsSolverOptions(MechanicsSolverOptions):
    def __init__(
            self,
            forceSummary: bool=False,
            cylindricalStress: bool=False,
            sphericalStress: bool=False,
            RhieChowCorrection: bool=True,
            strainTensor: str="EulerAlmansi"
        ):
        super().__init__(
            type="largeStrainTotLag",
            forceSummary=forceSummary,
            cylindricalStress=cylindricalStress,
            sphericalStress=sphericalStress,
            RhieChowCorrection=RhieChowCorrection
        )

        self.strainTensor = strainTensor


    @property
    def strainTensor(self):
        return self._strainTensor

    @strainTensor.setter
    def strainTensor(self, strainTensor) -> None:
        check_type("strainTensor", strainTensor, str)
        self._strainTensor = strainTensor
        self.__setitem__('strainTensor', self.strainTensor)


class LargeStrainUpdLagMechanicsSolverOptions(MechanicsSolverOptions):
    def __init__(
            self,
            forceSummary: bool=False,
            cylindricalStress: bool=False,
            sphericalStress: bool=False,
            RhieChowCorrection: bool=True,
            strainTensor: str="EulerAlmansi"
        ):
        super().__init__(
            type="largeStrainUpdLag",
            forceSummary=forceSummary,
            cylindricalStress=cylindricalStress,
            sphericalStress=sphericalStress,
            RhieChowCorrection=RhieChowCorrection
        )

        self.strainTensor = strainTensor


    @property
    def strainTensor(self):
        return self._strainTensor

    @strainTensor.setter
    def strainTensor(self, strainTensor) -> None:
        check_type("strainTensor", strainTensor, str)
        self._strainTensor = strainTensor
        self.__setitem__('strainTensor', self.strainTensor)


#==============================================================================*
# Neutronics solver

class NeutronicsSolverOptions(OpenFOAMDict):
    def __init__(
            self,
            type: str="fromLatestTime"
        ):
        super().__init__()

        self.type = type


    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _OFFBEAT_NEUTRONICS_SOLVER_TYPES)
        self._type = type


class DiffusionNeutronicsSolverOptions(NeutronicsSolverOptions):
    def __init__(
            self
        ):
        super().__init__(type="diffusion")


#==============================================================================*
# Element transport solver

class ElementTransportSolverOptions(OpenFOAMDict):
    def __init__(
            self,
            type: str="fromLatestTime"
        ):
        super().__init__()

        self.type = type


    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _OFFBEAT_ELEMENT_TRANSPORT_SOLVER_TYPES)
        self._type = type


class ByListElementTransportSolverOptions(ElementTransportSolverOptions):
    def __init__(
            self,
            solvers: list[str] | List[str]
        ):
        super().__init__(type="byList")

        self.solvers = solvers


    @property
    def solvers(self):
        return self._solvers

    @solvers.setter
    def solvers(self, solvers) -> None:
        check_type("solvers", solvers, (list, List))
        if (isinstance(solvers, List)):
            self._solvers = solvers
        else:
            self._solvers = List(solvers)
        self.__setitem__('solvers', self.solvers)


#==============================================================================*
# Burnup

class BurnupOptions(OpenFOAMDict):
    def __init__(
            self,
            type: str="none"
        ):
        super().__init__()
        self.type = type


    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _OFFBEAT_BURNUP_TYPES)
        self._type = type


class ConstantBurnupOptions(BurnupOptions):
    def __init__(
            self
        ):
        super().__init__(type="fromLatestTime")


class FromPowerBurnupOptions(BurnupOptions):
    def __init__(
            self,
            heatSourceName: str="Q",
        ):
        super().__init__(type="fromPower")

        self.heatSourceName = heatSourceName


    @property
    def heatSourceName(self):
        return self._heatSourceName

    @heatSourceName.setter
    def heatSourceName(self, heatSourceName) -> None:
        check_type("heatSourceName", heatSourceName, str)
        self._heatSourceName = heatSourceName
        self.__setitem__('heatSourceName', self.heatSourceName)


class LassmannBurnupOptions(BurnupOptions):
    def __init__(
            self,
            convergencePrecision: float=1e-2,
        ):
        super().__init__(type="Lassmann")

        self.convergencePrecision = convergencePrecision


    @property
    def convergencePrecision(self):
        return self._convergencePrecision

    @convergencePrecision.setter
    def convergencePrecision(self, convergencePrecision) -> None:
        check_type("convergencePrecision", convergencePrecision, (float, int))
        self._convergencePrecision = convergencePrecision
        self.__setitem__('convergencePrecision', self.convergencePrecision)


#==============================================================================*
# Fission Gas Release

class FgrOptions(OpenFOAMDict):
    def __init__(
            self,
            type: str="none"
        ):
        super().__init__()
        self.type = type


    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _OFFBEAT_FGR_TYPES)
        self._type = type


class SciantixFgrOptions(FgrOptions):
    def __init__(
            self,
            nFrequency: float=1,
            relax: float=1,
        ):
        super().__init__(type="SCIANTIX")

        self.nFrequency = nFrequency
        self.relax = relax


    @property
    def nFrequency(self):
        return self._nFrequency

    @nFrequency.setter
    def nFrequency(self, nFrequency) -> None:
        check_type("nFrequency", nFrequency, (float, int))
        self._nFrequency = nFrequency
        self.__setitem__('nFrequency', self.nFrequency)

    @property
    def relax(self):
        return self._relax

    @relax.setter
    def relax(self, relax) -> None:
        check_type("relax", relax, (float, int))
        self._relax = relax
        self.__setitem__('relax', self.relax)


#==============================================================================*
# Gap Gas

class GapGasOptions(OpenFOAMDict):
    def __init__(
            self,
            type: str
        ):
        super().__init__()
        self.type = type


    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _OFFBEAT_GAP_GAS_TYPES)
        self._type = type


class TrisoGapGasOptions(GapGasOptions):
    def __init__(
            self,
            bufferPorosity: float,
            COProductionModel: str,
            bufferOuterPatches: str,
            gapVolumeOffset: float=0,
            gasReserveVolume: float=0,
            gasReserveTemperature: float=290,
        ):
        super().__init__(type="TRISO")

        self.bufferPorosity = bufferPorosity
        self.COProductionModel = COProductionModel
        self.bufferOuterPatches = bufferOuterPatches
        self.gapVolumeOffset = gapVolumeOffset
        self.gasReserveVolume = gasReserveVolume
        self.gasReserveTemperature = gasReserveTemperature


    def __repr__(self, depth: int=0):
        self.__setitem__("bufferPorosity", self.bufferPorosity)
        self.__setitem__("COProductionModel", self.COProductionModel)
        self.__setitem__("bufferOuterPatches", self.bufferOuterPatches)

        return super().__repr__(depth)


    @property
    def bufferPorosity(self):
        return self._bufferPorosity

    @bufferPorosity.setter
    def bufferPorosity(self, bufferPorosity) -> None:
        check_type("bufferPorosity", bufferPorosity, (float, int))
        self._bufferPorosity = bufferPorosity

    @property
    def COProductionModel(self):
        return self._COProductionModel

    @COProductionModel.setter
    def COProductionModel(self, COProductionModel) -> None:
        check_type("COProductionModel", COProductionModel, str)
        self._COProductionModel = COProductionModel

    @property
    def bufferOuterPatches(self):
        return self._bufferOuterPatches

    @bufferOuterPatches.setter
    def bufferOuterPatches(self, bufferOuterPatches) -> None:
        check_type("bufferOuterPatches", bufferOuterPatches, str)
        self._bufferOuterPatches = bufferOuterPatches

    @property
    def gapVolumeOffset(self):
        return self._gapVolumeOffset

    @gapVolumeOffset.setter
    def gapVolumeOffset(self, gapVolumeOffset) -> None:
        check_type("gapVolumeOffset", gapVolumeOffset, (float, int))
        self._gapVolumeOffset = gapVolumeOffset
        self.__setitem__('gapVolumeOffset', self.gapVolumeOffset)

    @property
    def gasReserveVolume(self):
        return self._gasReserveVolume

    @gasReserveVolume.setter
    def gasReserveVolume(self, gasReserveVolume) -> None:
        check_type("gasReserveVolume", gasReserveVolume, (float, int))
        self._gasReserveVolume = gasReserveVolume
        self.__setitem__('gasReserveVolume', self.gasReserveVolume)

    @property
    def gasReserveTemperature(self):
        return self._gasReserveTemperature

    @gasReserveTemperature.setter
    def gasReserveTemperature(self, gasReserveTemperature) -> None:
        check_type("gasReserveTemperature", gasReserveTemperature, (float, int))
        self._gasReserveTemperature = gasReserveTemperature
        self.__setitem__('gasReserveTemperature', self.gasReserveTemperature)


class FrapconGapGasOptions(GapGasOptions):
    def __init__(
            self,
            gapVolumeOffset: float=0,
            gasReserveVolume: float=0,
            gasReserveTemperature: float=290,
            gapPatches: list | List=[],
            holePatches: list | List=[],
            topFuelPatches: list | List=[],
            bottomFuelPatches: list | List=[]
        ):
        super().__init__(type="FRAPCON")

        self.gapVolumeOffset = gapVolumeOffset
        self.gasReserveVolume = gasReserveVolume
        self.gasReserveTemperature = gasReserveTemperature
        self.gapPatches = gapPatches
        self.holePatches = holePatches
        self.topFuelPatches = topFuelPatches
        self.bottomFuelPatches = bottomFuelPatches


    def __repr__(self, depth = 0):
        self.__setitem__("gapPatches", self.gapPatches)
        self.__setitem__("holePatches", self.holePatches)
        self.__setitem__("topFuelPatches", self.topFuelPatches)
        self.__setitem__("bottomFuelPatches", self.bottomFuelPatches)

        return super().__repr__(depth)


    @property
    def gapVolumeOffset(self):
        return self._gapVolumeOffset

    @gapVolumeOffset.setter
    def gapVolumeOffset(self, gapVolumeOffset) -> None:
        check_type("gapVolumeOffset", gapVolumeOffset, (float, int))
        self._gapVolumeOffset = gapVolumeOffset
        self.__setitem__('gapVolumeOffset', self.gapVolumeOffset)

    @property
    def gasReserveVolume(self):
        return self._gasReserveVolume

    @gasReserveVolume.setter
    def gasReserveVolume(self, gasReserveVolume) -> None:
        check_type("gasReserveVolume", gasReserveVolume, (float, int))
        self._gasReserveVolume = gasReserveVolume
        self.__setitem__('gasReserveVolume', self.gasReserveVolume)

    @property
    def gasReserveTemperature(self):
        return self._gasReserveTemperature

    @gasReserveTemperature.setter
    def gasReserveTemperature(self, gasReserveTemperature) -> None:
        check_type("gasReserveTemperature", gasReserveTemperature, (float, int))
        self._gasReserveTemperature = gasReserveTemperature
        self.__setitem__('gasReserveTemperature', self.gasReserveTemperature)

    @property
    def gapPatches(self):
        return self._gapPatches

    @gapPatches.setter
    def gapPatches(self, gapPatches) -> None:
        check_type("gapPatches", gapPatches, (list, List, np.ndarray))
        if (isinstance(gapPatches, (list, np.ndarray))):
            self._gapPatches = List(gapPatches)
        else:
            self._gapPatches = gapPatches

    @property
    def holePatches(self):
        return self._holePatches

    @holePatches.setter
    def holePatches(self, holePatches) -> None:
        check_type("holePatches", holePatches, (list, List, np.ndarray), none_ok=True)
        if (isinstance(holePatches, (list, np.ndarray))):
            self._holePatches = List(holePatches)
        else:
            self._holePatches = holePatches

    @property
    def topFuelPatches(self):
        return self._topFuelPatches

    @topFuelPatches.setter
    def topFuelPatches(self, topFuelPatches) -> None:
        check_type("topFuelPatches", topFuelPatches, (list, List, np.ndarray), none_ok=True)
        if (isinstance(topFuelPatches, (list, np.ndarray))):
            self._topFuelPatches = List(topFuelPatches)
        else:
            self._topFuelPatches = topFuelPatches

    @property
    def bottomFuelPatches(self):
        return self._bottomFuelPatches

    @bottomFuelPatches.setter
    def bottomFuelPatches(self, bottomFuelPatches) -> None:
        check_type("bottomFuelPatches", bottomFuelPatches, (list, List, np.ndarray), none_ok=True)
        if (isinstance(bottomFuelPatches, (list, np.ndarray))):
            self._bottomFuelPatches = List(bottomFuelPatches)
        else:
            self._bottomFuelPatches = bottomFuelPatches


#==============================================================================*
# Rheology

class RheologyOptions(OpenFOAMDict):
    def __init__(
            self,
            type: str
        ):
        super().__init__()
        self.type = type


    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _OFFBEAT_RHEOLOGY_TYPES)
        self._type = type


class ByMaterialRheologyOptions(RheologyOptions):
    def __init__(
            self,
            thermalExpansion: bool=True,
            modifiedPlaneStrain: bool=True,
            springModulus: float=3500,
            coolantPressureList: Table=None,
            outOfBounds: str="clamp",
            planeStress: bool=False
        ):
        super().__init__(type="byMaterial")

        self.thermalExpansion = thermalExpansion
        self.modifiedPlaneStrain = modifiedPlaneStrain
        self.springModulus = springModulus
        self.coolantPressureList = coolantPressureList
        self.outOfBounds = outOfBounds
        self.planeStress = planeStress

    def __repr__(self, depth: int=0):
        self.coolantPressureList.type = ""
        self.__setitem__("coolantPressureList", OpenFOAMDict({
            "values": self.coolantPressureList,
            "outOfBounds": self.outOfBounds
        }))

        return super().__repr__(depth)


    @property
    def thermalExpansion(self):
        return self._thermalExpansion

    @thermalExpansion.setter
    def thermalExpansion(self, thermalExpansion) -> None:
        check_type("thermalExpansion", thermalExpansion, bool)
        self._thermalExpansion = thermalExpansion
        self.__setitem__('thermalExpansion', self.thermalExpansion)

    @property
    def modifiedPlaneStrain(self):
        return self._modifiedPlaneStrain

    @modifiedPlaneStrain.setter
    def modifiedPlaneStrain(self, modifiedPlaneStrain) -> None:
        check_type("modifiedPlaneStrain", modifiedPlaneStrain, bool)
        self._modifiedPlaneStrain = modifiedPlaneStrain
        self.__setitem__('modifiedPlaneStrain', self.modifiedPlaneStrain)

    @property
    def springModulus(self):
        return self._springModulus

    @springModulus.setter
    def springModulus(self, springModulus) -> None:
        check_type("springModulus", springModulus, (float, int))
        self._springModulus = springModulus
        self.__setitem__('springModulus', self.springModulus)

    @property
    def coolantPressureList(self):
        return self._coolantPressureList

    @coolantPressureList.setter
    def coolantPressureList(self, coolantPressureList) -> None:
        check_type("coolantPressureList", coolantPressureList, Table, none_ok=True)
        self._coolantPressureList = coolantPressureList

    @property
    def outOfBounds(self):
        return self._outOfBounds

    @outOfBounds.setter
    def outOfBounds(self, outOfBounds) -> None:
        check_type("outOfBounds", outOfBounds, str)
        self._outOfBounds = outOfBounds

    @property
    def planeStress(self):
        return self._planeStress

    @planeStress.setter
    def planeStress(self, planeStress) -> None:
        check_type("planeStress", planeStress, bool)
        self._planeStress = planeStress
        self.__setitem__('planeStress', self.planeStress)


#==============================================================================*
# Stress Analysis

class StressAnalysis(OpenFOAMDict):
    """
    StressAnalysis object

    Parameters
    ----------
    nCorrectors : int
    maxOuterIter : int
    referencePairs : List | list | np.ndarray
    D : Vector
    T : float
    relD : float
    relT : float
    """
    def __init__(
            self,
            nCorrectors: int=20,
            maxOuterIter: int=100,
            referencePairs: List | list | np.ndarray=[],
            D: Vector=Vector(1e-5, 0, 1e-5),
            T: float=1e-5,
            neutronFlux: float=1e-5,
            relD: float=1e-5,
            relT: float=1e-5,
            relneutronFlux0: float=1e-5
        ):
        super().__init__()
        self.nCorrectors = nCorrectors
        self.maxOuterIter = maxOuterIter
        self.referencePairs = referencePairs
        self.D = D
        self.T = T
        self.neutronFlux = neutronFlux
        self.relD = relD
        self.relT = relT
        self.relneutronFlux0 = relneutronFlux0

    @property
    def nCorrectors(self):
        return self._nCorrectors

    @nCorrectors.setter
    def nCorrectors(self, nCorrectors) -> None:
        check_type("nCorrectors", nCorrectors, int)
        self._nCorrectors = nCorrectors
        self.__setitem__('nCorrectors', self.nCorrectors)

    @property
    def maxOuterIter(self):
        return self._maxOuterIter

    @maxOuterIter.setter
    def maxOuterIter(self, maxOuterIter) -> None:
        check_type("maxOuterIter", maxOuterIter, int)
        self._maxOuterIter = maxOuterIter
        self.__setitem__('maxOuterIter', self.maxOuterIter)

    @property
    def referencePairs(self):
        return self._referencePairs

    @referencePairs.setter
    def referencePairs(self, referencePairs) -> None:
        check_type("referencePairs", referencePairs, (List, list, np.ndarray))
        if (isinstance(referencePairs, List)):
            self._referencePairs = referencePairs
        else:
            self._referencePairs = List(referencePairs)
        self.__setitem__('referencePairs', self.referencePairs)

    @property
    def D(self):
        return self._D

    @D.setter
    def D(self, D) -> None:
        check_type("D", D, Vector)
        self._D = D
        self.__setitem__('D', self.D)

    @property
    def T(self):
        return self._T

    @T.setter
    def T(self, T) -> None:
        check_type("T", T, (float, int))
        self._T = T
        self.__setitem__('T', self.T)

    @property
    def neutronFlux0(self):
        return self._neutronFlux0

    @neutronFlux0.setter
    def neutronFlux0(self, neutronFlux0) -> None:
        check_type("neutronFlux0", neutronFlux0, (float, int))
        self._neutronFlux0 = neutronFlux0
        self.__setitem__('neutronFlux0', self.neutronFlux0)

    @property
    def relD(self):
        return self._relD

    @relD.setter
    def relD(self, relD) -> None:
        check_type("relD", relD, (float, int))
        self._relD = relD
        self.__setitem__('relD', self.relD)

    @property
    def relT(self):
        return self._relT

    @relT.setter
    def relT(self, relT) -> None:
        check_type("relT", relT, (float, int))
        self._relT = relT
        self.__setitem__('relT', self.relT)

    @property
    def relneutronFlux0(self):
        return self._relneutronFlux0

    @relneutronFlux0.setter
    def relneutronFlux0(self, relneutronFlux0) -> None:
        check_type("relneutronFlux0", relneutronFlux0, (float, int))
        self._relneutronFlux0 = relneutronFlux0
        self.__setitem__('relneutronFlux0', self.relneutronFlux0)


#==============================================================================*
# Heat Source

class HeatSourceOptions(OpenFOAMDict):
    """
    Base class for heat source options.

    """
    def __init__(
            self,
            type: str="fromLatestTime",
            materials: list | List=[],
            axialProfile: AxialProfile=None,
            azimuthalProfile: AzimuthalProfile=None,
            radialProfile: RadialProfile=None,
            externalPowerUpdate: bool=None
        ):
        super().__init__()
        self.type = type
        self.materials = materials
        self.axialProfile = axialProfile
        self.azimuthalProfile = azimuthalProfile
        self.radialProfile = radialProfile
        self.externalPowerUpdate = externalPowerUpdate

    def __repr__(self, depth = 0):
        self.__setitem__('materials', self.materials)

        if (self.axialProfile is not None):
            self.__setitem__('axialProfile', self.axialProfile)
        if (self.azimuthalProfile is not None):
            self.__setitem__('azimuthalProfile', self.azimuthalProfile)
        if (self.radialProfile is not None):
            self.__setitem__('radialProfile', self.radialProfile)
        if (self.externalPowerUpdate is not None):
            self.__setitem__('externalPowerUpdate', self.externalPowerUpdate)

        return super().__repr__(depth)

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _OFFBEAT_HEAT_SOURCE_TYPES)
        self._type = type

    @property
    def materials(self):
        return self._materials

    @materials.setter
    def materials(self, materials) -> None:
        check_type("materials", materials, (list, List, np.ndarray))
        if (isinstance(materials, (list, np.ndarray))):
            self._materials = List(materials)
        else:
            self._materials = materials

    @property
    def axialProfile(self):
        return self._axialProfile

    @axialProfile.setter
    def axialProfile(self, axialProfile) -> None:
        check_type("axialProfile", axialProfile, AxialProfile, none_ok=True)
        self._axialProfile = axialProfile

    @property
    def azimuthalProfile(self):
        return self._azimuthalProfile

    @azimuthalProfile.setter
    def azimuthalProfile(self, azimuthalProfile) -> None:
        check_type("azimuthalProfile", azimuthalProfile, AzimuthalProfile, none_ok=True)
        self._azimuthalProfile = azimuthalProfile

    @property
    def radialProfile(self):
        return self._radialProfile

    @radialProfile.setter
    def radialProfile(self, radialProfile) -> None:
        check_type("radialProfile", radialProfile, RadialProfile, none_ok=True)
        self._radialProfile = radialProfile

    @property
    def externalPowerUpdate(self):
        return self._externalPowerUpdate

    @externalPowerUpdate.setter
    def externalPowerUpdate(self, externalPowerUpdate) -> None:
        check_type("externalPowerUpdate", externalPowerUpdate, bool, none_ok=True)
        self._externalPowerUpdate = externalPowerUpdate


class TimeDependentLhgrHeatSourceOptions(HeatSourceOptions):
    def __init__(
            self,
            timePoints: list | List=[],
            lhgr: list | List=[],
            timeInterpolationMethod: str="linear",
            materials: list | List=[],
            axialProfile: AxialProfile=None,
            azimuthalProfile: AzimuthalProfile=None,
            radialProfile: RadialProfile =None,
            externalPowerUpdate: bool=None
        ):
        super().__init__(
            type="timeDependentLhgr",
            materials=materials,
            axialProfile=axialProfile,
            azimuthalProfile=azimuthalProfile,
            radialProfile=radialProfile,
            externalPowerUpdate=externalPowerUpdate
        )

        self.timePoints = timePoints
        self.lhgr = lhgr
        self.timeInterpolationMethod = timeInterpolationMethod


    def __repr__(self, depth: int=0):
        self.__setitem__('timePoints', self.timePoints)
        self.__setitem__('lhgr', self.lhgr)
        self.__setitem__('timeInterpolationMethod', self.timeInterpolationMethod)

        return super().__repr__(depth)


    @property
    def timePoints(self):
        return self._timePoints

    @timePoints.setter
    def timePoints(self, timePoints) -> None:
        check_type("timePoints", timePoints, (list, List, np.ndarray))
        if (isinstance(timePoints, (list, np.ndarray))):
            self._timePoints = List(timePoints)
        else:
            self._timePoints = timePoints

    @property
    def lhgr(self):
        return self._lhgr

    @lhgr.setter
    def lhgr(self, lhgr) -> None:
        check_type("lhgr", lhgr, (list, List, np.ndarray))
        if (isinstance(lhgr, (list, np.ndarray))):
            self._lhgr = List(lhgr)
        else:
            self._lhgr = lhgr

    @property
    def timeInterpolationMethod(self):
        return self._timeInterpolationMethod

    @timeInterpolationMethod.setter
    def timeInterpolationMethod(self, timeInterpolationMethod) -> None:
        check_type("timeInterpolationMethod", timeInterpolationMethod, str)
        check_value("timeInterpolationMethod", timeInterpolationMethod, _INTERPOLATION_METHOD_TYPES)
        self._timeInterpolationMethod = timeInterpolationMethod


class TimeDependentVhgrHeatSourceOptions(HeatSourceOptions):
    def __init__(
            self,
            timePoints: list | List=[],
            vhgr: list | List=[],
            timeInterpolationMethod: str="linear",
            materials: list | List=[],
            axialProfile: AxialProfile=None,
            azimuthalProfile: AzimuthalProfile=None,
            radialProfile: RadialProfile =None,
            externalPowerUpdate: bool=None
        ):
        super().__init__(
            type="timeDependentVhgr",
            materials=materials,
            axialProfile=axialProfile,
            azimuthalProfile=azimuthalProfile,
            radialProfile=radialProfile,
            externalPowerUpdate=externalPowerUpdate
        )

        self.timePoints = timePoints
        self.vhgr = vhgr
        self.timeInterpolationMethod = timeInterpolationMethod


    def __repr__(self, depth: int=0):
        self.__setitem__('timePoints', self.timePoints)
        self.__setitem__('vhgr', self.vhgr)
        self.__setitem__('timeInterpolationMethod', self.timeInterpolationMethod)

        return super().__repr__(depth)


    @property
    def timePoints(self):
        return self._timePoints

    @timePoints.setter
    def timePoints(self, timePoints) -> None:
        check_type("timePoints", timePoints, (list, List, np.ndarray))
        if (isinstance(timePoints, (list, np.ndarray))):
            self._timePoints = List(timePoints)
        else:
            self._timePoints = timePoints

    @property
    def vhgr(self):
        return self._vhgr

    @vhgr.setter
    def vhgr(self, vhgr) -> None:
        check_type("vhgr", vhgr, (list, List, np.ndarray))
        if (isinstance(vhgr, (list, np.ndarray))):
            self._vhgr = List(vhgr)
        else:
            self._vhgr = vhgr

    @property
    def timeInterpolationMethod(self):
        return self._timeInterpolationMethod

    @timeInterpolationMethod.setter
    def timeInterpolationMethod(self, timeInterpolationMethod) -> None:
        check_type("timeInterpolationMethod", timeInterpolationMethod, str)
        check_value("timeInterpolationMethod", timeInterpolationMethod, _INTERPOLATION_METHOD_TYPES)
        self._timeInterpolationMethod = timeInterpolationMethod


#==============================================================================*
# Fast Flux

class FastFluxOptions(OpenFOAMDict):
    """
    Base class for fast flux options.

    """
    def __init__(
            self,
            type: str="fromLatestTime",
            materials: list | List=[],
            axialProfile: AxialProfile=None,
            radialProfile: RadialProfile=None
        ):
        super().__init__()
        self.type = type
        self.materials = materials
        self.axialProfile = axialProfile
        self.radialProfile = radialProfile

    def __repr__(self, depth: int=0):
        self.__setitem__('materials', self.materials)

        if (self.axialProfile is not None):
            self.__setitem__('axialProfile', self.axialProfile)
        if (self.radialProfile is not None):
            self.__setitem__('radialProfile', self.radialProfile)

        return super().__repr__(depth)


    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _OFFBEAT_FAST_FLUX_TYPES)
        self._type = type

    @property
    def materials(self):
        return self._materials

    @materials.setter
    def materials(self, materials) -> None:
        check_type("materials", materials, (list, List, np.ndarray))
        if (isinstance(materials, (list, np.ndarray))):
            self._materials = List(materials)
        else:
            self._materials = materials

    @property
    def axialProfile(self):
        return self._axialProfile

    @axialProfile.setter
    def axialProfile(self, axialProfile) -> None:
        check_type("axialProfile", axialProfile, AxialProfile, none_ok=True)
        self._axialProfile = axialProfile

    @property
    def radialProfile(self):
        return self._radialProfile

    @radialProfile.setter
    def radialProfile(self, radialProfile) -> None:
        check_type("radialProfile", radialProfile, RadialProfile, none_ok=True)
        self._radialProfile = radialProfile


class TimeDependentFastFluxOptions(FastFluxOptions):
    def __init__(
            self,
            timePoints: list | List=[],
            fastFlux: list | List=[],
            timeInterpolationMethod: str="linear",
            materials: list | List=[],
            axialProfile: AxialProfile=None,
            radialProfile: RadialProfile=None
        ):
        super().__init__(
            type="timeDependentAxialProfile",
            materials=materials,
            axialProfile=axialProfile,
            radialProfile=radialProfile
        )
        self.timePoints = timePoints
        self.fastFlux = fastFlux
        self.timeInterpolationMethod = timeInterpolationMethod


    def __repr__(self, depth: int=0):
        self.__setitem__('timePoints', self.timePoints)
        self.__setitem__('fastFlux', self.fastFlux)
        self.__setitem__('timeInterpolationMethod', self.timeInterpolationMethod)

        return super().__repr__(depth)


    @property
    def timePoints(self):
        return self._timePoints

    @timePoints.setter
    def timePoints(self, timePoints) -> None:
        check_type("timePoints", timePoints, (list, List, np.ndarray))
        if (isinstance(timePoints, (list, np.ndarray))):
            self._timePoints = List(timePoints)
        else:
            self._timePoints = timePoints

    @property
    def fastFlux(self):
        return self._fastFlux

    @fastFlux.setter
    def fastFlux(self, fastFlux) -> None:
        check_type("fastFlux", fastFlux, (list, List, np.ndarray))
        if (isinstance(fastFlux, (list, np.ndarray))):
            self._fastFlux = List(fastFlux)
        else:
            self._fastFlux = fastFlux

    @property
    def timeInterpolationMethod(self):
        return self._timeInterpolationMethod

    @timeInterpolationMethod.setter
    def timeInterpolationMethod(self, timeInterpolationMethod) -> None:
        check_type("timeInterpolationMethod", timeInterpolationMethod, str)
        check_value("timeInterpolationMethod", timeInterpolationMethod, _INTERPOLATION_METHOD_TYPES)
        self._timeInterpolationMethod = timeInterpolationMethod


#==============================================================================*
# Corrosion

class CorrosionOptions(OpenFOAMDict):
    def __init__(
            self,
            type: str="fromLatestTime"
        ):
        super().__init__()
        self.type = type


    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _OFFBEAT_CORROSION_TYPES)
        self._type = type


#==============================================================================*
# Slice mapper

class SliceMapperOptions(OpenFOAMDict):
    def __init__(
            self,
            type: str="none"
        ):
        super().__init__()
        self.type = type


    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _OFFBEAT_SLICE_MAPPER_TYPES)
        self._type = type


class ByMaterialSliceMapperOptions(SliceMapperOptions):
    def __init__(
            self
        ):
        super().__init__(type="byMaterial")


class AutoAxialSliceMapperOptions(SliceMapperOptions):
    def __init__(
            self,
            precision: float=1e-6
        ):
        super().__init__(type="autoAxialSlices")

        self.precision = precision


    @property
    def precision(self):
        return self._precision

    @precision.setter
    def precision(self, precision) -> None:
        check_type("precision", precision, (float, int))
        self._precision = precision
        self.__setitem__('precision', self.precision)


class ByPelletSliceMapperOptions(SliceMapperOptions):
    def __init__(
            self
        ):
        super().__init__(type="byPellets")


#==============================================================================*
# OFFBEAT Main Solver

class OffbeatSolver(Solver):
    """
    OFFBEAT solver parameter object.

    More info: https://gitlab.com/foam-for-nuclear/offbeat/-/tree/develop?ref_type=heads

    Parameters
    ----------
    region : str
        Name of the Offbeat region
    solver : str {"offbeat", "extendedThermoMechanics"}
        Offbeat solver name. If use in GeN-Foam, please use
        `extendedThermoMechanics`.
    thermalSolverOptions : ThermalSolverOptions

    mechanicsSolverOptions : MechanicsSolverOptions

    neutronicsSolver : NeutronicsSolverOptions

    elementTransportSolverOptions : ElementTransportSolverOptions

    materialProperties : str {"byZone"}
        Default "byZone"
    rheologyOptions : RheologyOptions

    heatSourceOptions : HeatSourceOptions

    burnupOptions : BurnupOptions

    fastFluxOptions : FastFluxOptions

    corrosionOptions : CorrosionOptions

    gapGasOptions : GapGasOptions

    fgrOptions : FgrOptions

    sliceMapperOptions : SliceMapperOptions

    globalOptions: GlobalOptions

    couplingOptions: ThermoMechanicsCouplingOptions

    removeBaffles : bool
        Default `False`
    timeFolder : TimeFolder
        Object representing the time folder.
    mesh : Mesh
        Mesh object
    isMeshDeformation : bool
        Flag to allow mesh deformation based on a displacement field
        (default `False`).
    displacementFieldName : str
        Name of the displacement field (default `disp`).
    isSetFvSolutionToDefault : bool
        Flag to set fvSolution with default parameters (default `True`).
    isSetFvSchemesToDefault : bool
        Flag to set fvSchemes with default parameters (default `True`).

    Attributes
    ----------
    couplingOptions : ThermoMechanicsCouplingOptions
    globalOptions : GlobalOptions
    thermalSolverOptions : ThermalSolverOptions
    mechanicsSolverOptions : MechanicsSolverOptions
    neutronicsSolverOptions : NeutronicsSolverOptions
    elementTransportSolverOptions : ElementTransportSolverOptions
    rheologyOptions: RheologyOptions
    fgrOptions : FgrOptions
    gapGasOptions : GapGasOptions
    heatSourceOptions : HeatSourceOptions
    fastFluxOptions : FastFluxOptions
    corrosionOptions : CorrosionOptions
    burnupOptions : BurnupOptions
    sliceMapperOptions : SliceMapperOptions
    materials : OpenFOAMListDict
    stressAnalysis : StressAnalysis
    """

    def __init__(
            self,
            solver: str,
            region: str="",
            thermalSolverOptions: ThermalSolverOptions=ThermalSolverOptions(),
            mechanicsSolverOptions: MechanicsSolverOptions=MechanicsSolverOptions(),
            neutronicsSolverOptions: NeutronicsSolverOptions=NeutronicsSolverOptions(),
            elementTransportSolverOptions: ElementTransportSolverOptions=ElementTransportSolverOptions(),
            materialProperties: str="byZone",
            rheologyOptions: RheologyOptions=ByMaterialRheologyOptions(),
            heatSourceOptions: HeatSourceOptions=HeatSourceOptions(),
            burnupOptions: BurnupOptions=BurnupOptions(),
            fastFluxOptions: FastFluxOptions=FastFluxOptions(),
            corrosionOptions: CorrosionOptions=CorrosionOptions(),
            gapGasOptions: str=FrapconGapGasOptions(),
            fgrOptions: FgrOptions=FgrOptions(),
            sliceMapperOptions: SliceMapperOptions=AutoAxialSliceMapperOptions(),
            globalOptions: GlobalOptions=GlobalOptions(),
            couplingOptions: ThermoMechanicsCouplingOptions=ThermoMechanicsCouplingOptions(),
            removeBaffles: bool=False,
            timeFolder: TimeFolder=None,
            mesh: Mesh=None,
            isMeshDeformation: bool=False,
            displacementFieldName: str="disp",
            isSetFvSolutionToDefault: bool=True,
            isSetFvSchemesToDefault: bool=True
        ):
        super().__init__(region, solver, removeBaffles, timeFolder, mesh, isMeshDeformation, displacementFieldName)

        self.materialProperties = materialProperties

        self.couplingOptions: ThermoMechanicsCouplingOptions = couplingOptions
        self.globalOptions: GlobalOptions = globalOptions
        self.thermalSolverOptions: ThermalSolverOptions = thermalSolverOptions
        self.mechanicsSolverOptions: MechanicsSolverOptions = mechanicsSolverOptions
        self.neutronicsSolverOptions: NeutronicsSolverOptions = neutronicsSolverOptions
        self.elementTransportSolverOptions: ElementTransportSolverOptions = elementTransportSolverOptions
        self.rheologyOptions: RheologyOptions = rheologyOptions
        self.fgrOptions: FgrOptions = fgrOptions
        self.gapGasOptions: GapGasOptions = gapGasOptions
        self.heatSourceOptions: HeatSourceOptions = heatSourceOptions
        self.fastFluxOptions: FastFluxOptions = fastFluxOptions
        self.corrosionOptions: CorrosionOptions = corrosionOptions
        self.burnupOptions: BurnupOptions = burnupOptions
        self.sliceMapperOptions: SliceMapperOptions = sliceMapperOptions

        self.materials = OpenFOAMListDict(BaseThermomechanicalMaterial, "materials")

        self.stressAnalysis: StressAnalysis = StressAnalysis()

        if (isSetFvSolutionToDefault):
            self.set_fvSolution_default()
        if (isSetFvSchemesToDefault):
            self.set_fvSchemes_default()


    @property
    def solver(self):
        return self._solver

    @solver.setter
    def solver(self, solver) -> None:
        check_type("solver", solver, str)
        check_value("solver", solver, _OFFBEAT_SOLVER_TYPES)
        self._solver = solver

    @property
    def thermalSolverOptions(self):
        return self._thermalSolverOptions

    @thermalSolverOptions.setter
    def thermalSolverOptions(self, thermalSolverOptions) -> None:
        check_type("thermalSolverOptions", thermalSolverOptions, ThermalSolverOptions)
        self._thermalSolverOptions = thermalSolverOptions

    @property
    def mechanicsSolverOptions(self):
        return self._mechanicsSolverOptions

    @mechanicsSolverOptions.setter
    def mechanicsSolverOptions(self, mechanicsSolverOptions) -> None:
        check_type("mechanicsSolverOptions", mechanicsSolverOptions, MechanicsSolverOptions)
        self._mechanicsSolverOptions = mechanicsSolverOptions

    @property
    def neutronicsSolverOptions(self):
        return self._neutronicsSolverOptions

    @neutronicsSolverOptions.setter
    def neutronicsSolverOptions(self, neutronicsSolverOptions) -> None:
        check_type("neutronicsSolverOptions", neutronicsSolverOptions, NeutronicsSolverOptions)
        self._neutronicsSolverOptions = neutronicsSolverOptions

    @property
    def elementTransportSolverOptions(self):
        return self._elementTransportSolverOptions

    @elementTransportSolverOptions.setter
    def elementTransportSolverOptions(self, elementTransportSolverOptions) -> None:
        check_type("elementTransportSolverOptions", elementTransportSolverOptions, ElementTransportSolverOptions)
        self._elementTransportSolverOptions = elementTransportSolverOptions

    @property
    def materialProperties(self):
        return self._materialProperties

    @materialProperties.setter
    def materialProperties(self, materialProperties) -> None:
        check_type("materialProperties", materialProperties, str)
        check_value("materialProperties", materialProperties, _OFFBEAT_MATERIAL_PROPERTIES_TYPES)
        self._materialProperties = materialProperties

    @property
    def rheologyOptions(self):
        return self._rheologyOptions

    @rheologyOptions.setter
    def rheologyOptions(self, rheologyOptions) -> None:
        check_type("rheologyOptions", rheologyOptions, RheologyOptions)
        self._rheologyOptions = rheologyOptions

    @property
    def burnupOptions(self):
        return self._burnupOptions

    @burnupOptions.setter
    def burnupOptions(self, burnupOptions) -> None:
        check_type("burnupOptions", burnupOptions, BurnupOptions)
        self._burnupOptions = burnupOptions

    @property
    def corrosionOptions(self):
        return self._corrosionOptions

    @corrosionOptions.setter
    def corrosionOptions(self, corrosionOptions) -> None:
        check_type("corrosionOptions", corrosionOptions, CorrosionOptions)
        self._corrosionOptions = corrosionOptions

    @property
    def fgrOptions(self):
        return self._fgrOptions

    @fgrOptions.setter
    def fgrOptions(self, fgrOptions) -> None:
        check_type("fgrOptions", fgrOptions, FgrOptions)
        self._fgrOptions = fgrOptions

    @property
    def sliceMapperOptions(self):
        return self._sliceMapperOptions

    @sliceMapperOptions.setter
    def sliceMapperOptions(self, sliceMapperOptions) -> None:
        check_type("sliceMapperOptions", sliceMapperOptions, SliceMapperOptions)
        self._sliceMapperOptions = sliceMapperOptions

    @property
    def stressAnalysis(self):
        return self._stressAnalysis

    @stressAnalysis.setter
    def stressAnalysis(self, stressAnalysis) -> None:
        check_type("stressAnalysis", stressAnalysis, StressAnalysis)
        self._stressAnalysis = stressAnalysis

    @property
    def gapGasOptions(self):
        return self._gapGasOptions

    @gapGasOptions.setter
    def gapGasOptions(self, gapGasOptions) -> None:
        check_type("gapGasOptions", gapGasOptions, GapGasOptions)
        self._gapGasOptions = gapGasOptions

    @property
    def heatSourceOptions(self):
        return self._heatSourceOptions

    @heatSourceOptions.setter
    def heatSourceOptions(self, heatSourceOptions) -> None:
        check_type("heatSourceOptions", heatSourceOptions, HeatSourceOptions)
        self._heatSourceOptions = heatSourceOptions

    @property
    def fastFluxOptions(self):
        return self._fastFluxOptions

    @fastFluxOptions.setter
    def fastFluxOptions(self, fastFluxOptions) -> None:
        check_type("fastFluxOptions", fastFluxOptions, FastFluxOptions)
        self._fastFluxOptions = fastFluxOptions

    @property
    def is_offbeat_solver(self) -> bool:
        return(self.solver == "offbeat")

    @property
    def is_extended_thermomechanics_solver(self) -> bool:
        return(self.solver == "extendedThermoMechanics")


    def add_material(self, material: BaseThermomechanicalMaterial):
        check_type("material", material, BaseThermomechanicalMaterial)
        self.materials.append(material)


    def set_fvSolution_default(self):
        self.fvSolution = fvSolution()

        self.fvSolution.solvers["porosity"] = fvSolutionSolver(
            solver="PBiCG",
            smoother="GaussSeidel",
            preconditioner="DILU",
            tolerance=1e-10,
            relTol=1e-2,
        )
        self.fvSolution.solvers["D"] = fvSolutionSolver(
            solver="PCG",
            preconditioner="FDIC",
            tolerance=1e-10,
            relTol=0.1,
        )
        self.fvSolution.solvers["T"] = copy(self.fvSolution.solvers["D"])
        self.fvSolution.solvers["neutronFlux0"] = copy(self.fvSolution.solvers["D"])

        self.fvSolution.solvers["fuelDisp"] = fvSolutionSolver(
            solver="PBiCG",
            preconditioner="DILU",
            tolerance=1e-5,
            relTol=1e-2,
        )
        self.fvSolution.solvers["CRDisp"] = copy(self.fvSolution.solvers["fuelDisp"])

        self.stressAnalysis.nCorrectors = 20
        self.stressAnalysis.maxOuterIter = 100
        self.stressAnalysis.referencePairs = []
        self.stressAnalysis.D = Vector(1e-5, 1, 1e-5)
        self.stressAnalysis.T = 1e-5
        self.stressAnalysis.relD = 1e-5
        self.stressAnalysis.relT = 1e-5

        self.add_relaxation_on_field('D', 0.7)


    def set_fvSchemes_default(self):
        self.fvSchemes = fvSchemes()

        self.fvSchemes.d2dt2Schemes['default'] = 'Euler'
        self.fvSchemes.ddtSchemes['default'] = 'Euler'

        self.fvSchemes.gradSchemes['default'] = 'leastSquaresImplicitContact'

        self.fvSchemes.divSchemes['default'] = 'Gauss linear'
        for scheme in [
            'div(flux,P)', 'div(fuelDisp)', 'div(CRDisp)'
        ]:
            self.fvSchemes.divSchemes[scheme] = 'Gauss upwind'

        self.fvSchemes.laplacianSchemes['default'] = 'Gauss linear corrected'

        self.fvSchemes.interpolationSchemes['default'] = 'linear'

        self.fvSchemes.snGradSchemes['default'] = 'uncorrected'

        self.fvSchemes.fluxRequired['default'] = False
        self.fvSchemes.fluxRequired['D'] = None
        self.fvSchemes.fluxRequired['T'] = None


    def get_required_fields(self) -> list[str]:
        """
        Return a list of required fields depending on the solver selected.
        They are labeled "MUST_READ"
        """
        solverFields = ["T"]
        if (self.mechanicsSolverOptions.type != "fromLatestTime"):
            solverFields += ["D"]
        if (self.is_offbeat_solver):
            solverFields += []
        if (self.is_extended_thermomechanics_solver):
            if (self.couplingOptions.correctDispForNeutro):
                solverFields += ["fuelDisp", "CRDisp"]

        return(solverFields)


    def get_default_fields(self) -> list[str]:
        """
        Return a list of default fields that can be used as input depending on
        the solver selected. They are labeled "READ_IF_PRESENT"
        """
        common = []
        solverSpecific = []

        if (self.is_offbeat_solver):
            solverSpecific = []
        elif (self.is_extended_thermomechanics_solver):
            solverSpecific = []

        return(self.get_required_fields() + solverSpecific + common)


    def export_properties_to_openfoam(self):
        # Dict file name
        filename = ""
        if (self.is_offbeat_solver):
            filename = "solverDict"
        elif (self.is_extended_thermomechanics_solver):
            filename = "thermoMechanicalProperties"

        # Create the dict
        text = openfoamHeader
        text += openfoamFileHeader(filename)

        text += addParameter('thermalSolver', self.thermalSolverOptions.type)
        text += addParameter('mechanicsSolver', self.mechanicsSolverOptions.type)
        if (self.is_offbeat_solver):
            text += addParameter('neutronicsSolver', self.neutronicsSolverOptions.type)
            text += addParameter('elementTransport', self.elementTransportSolverOptions.type, isAddExtraLine=True)

        text += addParameter('materialProperties', self.materialProperties)
        text += addParameter('rheology', self.rheologyOptions.type, isAddExtraLine=True)

        text += addParameter('heatSource', self.heatSourceOptions.type)
        text += addParameter('burnup', self.burnupOptions.type)
        text += addParameter('fastFlux', self.fastFluxOptions.type)
        text += addParameter('gapGas', self.gapGasOptions.type)
        text += addParameter('fgr', self.fgrOptions.type)
        text += addParameter('sliceMapper', self.sliceMapperOptions.type)
        text += addParameter('corrosion', self.corrosionOptions.type, isAddExtraLine=True)


        if (self.is_extended_thermomechanics_solver):
            text += f"couplingOptions{self.couplingOptions!r}\n"

        text += f"globalOptions{self.globalOptions!r}\n"
        text += f"thermalSolverOptions{self.thermalSolverOptions!r}\n"
        text += f"mechanicsSolverOptions{self.mechanicsSolverOptions!r}\n"
        if (self.is_offbeat_solver):
            text += f"rheologyOptions{self.rheologyOptions!r}\n"
            text += f"fgrOptions{self.fgrOptions!r}\n"
            text += f"gapGasOptions{self.gapGasOptions!r}\n"
            text += f"heatSourceOptions{self.heatSourceOptions!r}\n"
            text += f"fastFluxOptions{self.fastFluxOptions!r}\n"

            if (not self.neutronicsSolverOptions.is_empty):
                text += f"neutronicsSolverOption{self.neutronicsSolverOptions!r}\n"
            if (not self.elementTransportSolverOptions.is_empty):
                text += f"elementTransportOptions{self.elementTransportSolverOptions!r}\n"
            if (not self.corrosionOptions.is_empty):
                text += f"corrosionOptions{self.corrosionOptions!r}\n"
            if (not self.burnupOptions.is_empty):
                text += f"burnupOptions{self.burnupOptions!r}\n"
            if (not self.sliceMapperOptions.is_empty):
                text += f"mapperOptions{self.sliceMapperOptions!r}\n"

        text += f"{self.materials!r}"

        text += openfoamFooterLine

        # Write the file at the end
        with open(f"constant/{self.region}/{filename}", 'w') as f:
            f.write(text)


    def export_to_openfoam(self):
        self.create_folders()

        self.export_properties_to_openfoam()

        self.fvSolution.extraDict['stressAnalysis'] = self.stressAnalysis

        return super().export_to_openfoam()
