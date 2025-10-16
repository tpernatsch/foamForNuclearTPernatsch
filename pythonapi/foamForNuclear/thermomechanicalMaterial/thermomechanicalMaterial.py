from foamForNuclear.checkvalue import check_positive, check_type, check_value
from foamForNuclear.common import OpenFOAMDict, Table
from foamForNuclear.field import ReducedDimension


_OFFBEAT_MATERIAL_TYPES = {
    "constant",
    "UO2", "UPuO2",
    "zircaloy", "steel1515Ti", "HastelloyN", "inconel600", "molybdenum", "PyC",
    "SiC"
}


class BaseThermomechanicalMaterial(OpenFOAMDict):
    def __init__(
            self,
            name: str,
            material: str
        ):
        super().__init__(name=name)

        self.material = material


    @property
    def material(self):
        return self._material

    @material.setter
    def material(self, material) -> None:
        check_type("material", material, str)
        check_value("material", material, _OFFBEAT_MATERIAL_TYPES)
        self._material = material
        self.__setitem__('material', self.material)


class ThermomechanicalPropertyModel:
    def __init__(
            self,
            name: str
        ):
        self.name = name


    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name) -> None:
        check_type("name", name, str)
        self._name = name


class ConductivityModel(ThermomechanicalPropertyModel):
    def __init__(
            self,
            name: str
        ):
        super().__init__(name=name)

    def setitem(self, parentMaterial: OpenFOAMDict):
        parentMaterial.__setitem__("conductivityModel", self.name)


class DensityModel(ThermomechanicalPropertyModel):
    def __init__(
            self,
            name: str
        ):
        super().__init__(name=name)

    def setitem(self, parentMaterial: OpenFOAMDict):
        parentMaterial.__setitem__("densityModel", self.name)


class EmissivityModel(ThermomechanicalPropertyModel):
    def __init__(
            self,
            name: str
        ):
        super().__init__(name=name)

    def setitem(self, parentMaterial: OpenFOAMDict):
        parentMaterial.__setitem__("emissivityModel", self.name)


class HeatCapacityModel(ThermomechanicalPropertyModel):
    def __init__(
            self,
            name: str
        ):
        super().__init__(name=name)

    def setitem(self, parentMaterial: OpenFOAMDict):
        parentMaterial.__setitem__("heatTransferModel", self.name)


class PoissonRatioModel(ThermomechanicalPropertyModel):
    def __init__(
            self,
            name: str
        ):
        super().__init__(name=name)

    def setitem(self, parentMaterial: OpenFOAMDict):
        parentMaterial.__setitem__("PoissonRatioModel", self.name)


class ThermalExpansionModel(ThermomechanicalPropertyModel):
    def __init__(
            self,
            name: str
        ):
        super().__init__(name=name)

    def setitem(self, parentMaterial: OpenFOAMDict):
        parentMaterial.__setitem__("thermalExpansionModel", self.name)


class YoungModulusModel(ThermomechanicalPropertyModel):
    def __init__(
            self,
            name: str
        ):
        super().__init__(name=name)

    def setitem(self, parentMaterial: OpenFOAMDict):
        parentMaterial.__setitem__("YoungModulusModel", self.name)


class DensificationModel(ThermomechanicalPropertyModel):
    def __init__(
            self,
            name: str
        ):
        super().__init__(name=name)

    def setitem(self, parentMaterial: OpenFOAMDict):
        parentMaterial.__setitem__("densificationModel", self.name)


class SwellingModel(ThermomechanicalPropertyModel):
    def __init__(
            self,
            name: str
        ):
        super().__init__(name=name)

    def setitem(self, parentMaterial: OpenFOAMDict):
        parentMaterial.__setitem__("swellingModel", self.name)


class PhaseTransitionModel(ThermomechanicalPropertyModel):
    def __init__(
            self,
            name: str
        ):
        super().__init__(name=name)

    def setitem(self, parentMaterial: OpenFOAMDict):
        parentMaterial.__setitem__("phaseTransitionModel", self.name)


class RelocationModel(ThermomechanicalPropertyModel):
    def __init__(
            self,
            name: str
        ):
        super().__init__(name=name)

    def setitem(self, parentMaterial: OpenFOAMDict):
        parentMaterial.__setitem__("relocationModel", self.name)


class FailureModel(ThermomechanicalPropertyModel):
    def __init__(
            self,
            name: str="none"
        ):
        super().__init__(name=name)

    def setitem(self, parentMaterial: OpenFOAMDict):
        parentMaterial.__setitem__("failureModel", self.name)


class PoreVelocityModel(ThermomechanicalPropertyModel):
    def __init__(
            self,
            name: str
        ):
        super().__init__(name=name)

    def setitem(self, parentMaterial: OpenFOAMDict):
        parentMaterial.__setitem__("poreVelocityModel", self.name)


class YieldStressModel():
    def __init__(self, name: str="hardening"):
        self.name = name

    def setitem(self, parentMaterial: OpenFOAMDict):
        parentMaterial.__setitem__("yieldStressModel", self.name)


class HardningYieldStressModel(YieldStressModel):
    def __init__(
            self,
            plasticStrainVsYieldStress: Table
        ):
        super().__init__(name="hardening")

        self.plasticStrainVsYieldStress = plasticStrainVsYieldStress


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("plasticStrainVsYieldStress", self.plasticStrainVsYieldStress)


class CreepModel():
    def __init__(self, name: str):
        self.name = name

    def setitem(self, parentMaterial: OpenFOAMDict):
        parentMaterial.__setitem__("creepModel", self.name)


class LimbackCreepModel(CreepModel):
    def __init__(
            self,
            relax: float=1.0
        ):
        super().__init__(name="LimbackCreepModel")

        self.relax = relax


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("relax", self.relax)


class RheologyConstitutiveLaw(ThermomechanicalPropertyModel):
    def __init__(
            self,
            name: str,
            yieldStressModel: YieldStressModel=None,
            creepModel: CreepModel=None
        ):
        super().__init__(name)

        self.rheologyModelOptions = OpenFOAMDict({})

        self.yieldStressModel: YieldStressModel = yieldStressModel
        self.creepModel: CreepModel = creepModel


    def setitem(self, parentMaterial: OpenFOAMDict):
        parentMaterial.__setitem__("rheologyModel", self.name)

        if (self.yieldStressModel is not None):
            self.yieldStressModel.setitem(self.rheologyModelOptions)

        if (self.creepModel is not None):
            self.creepModel.setitem(self.rheologyModelOptions)

        if (not self.rheologyModelOptions.is_empty):
            parentMaterial.__setitem__("rheologyModelOptions", self.rheologyModelOptions)

    @property
    def yieldStressModel(self):
        return self._yieldStressModel

    @yieldStressModel.setter
    def yieldStressModel(self, yieldStressModel) -> None:
        check_type("yieldStressModel", yieldStressModel, YieldStressModel, none_ok=True)
        self._yieldStressModel = yieldStressModel

    @property
    def creepModel(self):
        return self._creepModel

    @creepModel.setter
    def creepModel(self, creepModel) -> None:
        check_type("creepModel", creepModel, CreepModel, none_ok=True)
        self._creepModel = creepModel


class ElasticityRheologyModel(RheologyConstitutiveLaw):
    def __init__(
            self,
            yieldStressModel: YieldStressModel=None,
            creepModel: CreepModel=None
        ):
        super().__init__(
            name="elasticity",
            yieldStressModel=yieldStressModel,
            creepModel=creepModel
        )


class MisesPlasticityRheologyModel(RheologyConstitutiveLaw):
    def __init__(
            self,
            yieldStressModel: YieldStressModel=None,
            creepModel: CreepModel=None
        ):
        super().__init__(
            name="misesPlasticity",
            yieldStressModel=yieldStressModel,
            creepModel=creepModel
        )


class MisesPlasticCreepRheologyModel(RheologyConstitutiveLaw):
    def __init__(
            self,
            yieldStressModel: YieldStressModel=None,
            creepModel: CreepModel=None
        ):
        super().__init__(
            name="misesPlasticCreep",
            yieldStressModel=yieldStressModel,
            creepModel=creepModel
        )


class MaterialModel(BaseThermomechanicalMaterial):
    def __init__(
            self,
            name: str,
            material: str,
            Tref: float,
            densityModel: DensityModel,
            heatCapacityModel: HeatCapacityModel,
            conductivityModel: ConductivityModel,
            emissivityModel: EmissivityModel,
            poissonRatioModel: PoissonRatioModel,
            thermalExpansionModel: ThermalExpansionModel,
            youngModulusModel: YoungModulusModel,
            densificationModel: DensificationModel=None,
            swellingModel: SwellingModel=None,
            phaseTransitionModel: PhaseTransitionModel=None,
            relocationModel: RelocationModel=None,
            failureModel: FailureModel=None,
            poreVelocityModel: PoreVelocityModel=None,
            rheologyModel: RheologyConstitutiveLaw=None
        ):
        super().__init__(name=name, material=material)

        self.Tref = Tref
        self.densityModel = densityModel
        self.heatCapacityModel = heatCapacityModel
        self.conductivityModel = conductivityModel
        self.emissivityModel = emissivityModel
        self.poissonRatioModel = poissonRatioModel
        self.thermalExpansionModel = thermalExpansionModel
        self.youngModulusModel = youngModulusModel

        self.densificationModel = densificationModel
        self.swellingModel = swellingModel
        self.phaseTransitionModel = phaseTransitionModel
        self.relocationModel = relocationModel
        self.failureModel = failureModel
        self.poreVelocityModel = poreVelocityModel

        self.rheologyModel = rheologyModel


    def __repr__(self, depth: int=0):
        self.__setitem__('Tref', ConstantMechanicalLaw("Tref", ReducedDimension('T'), self.Tref))

        self.densityModel.setitem(self)
        self.heatCapacityModel.setitem(self)
        self.conductivityModel.setitem(self)
        self.emissivityModel.setitem(self)
        self.poissonRatioModel.setitem(self)
        self.thermalExpansionModel.setitem(self)
        self.youngModulusModel.setitem(self)

        if (self.densificationModel is not None):
            self.densificationModel.setitem(self)
        if (self.swellingModel is not None):
            self.swellingModel.setitem(self)
        if (self.phaseTransitionModel is not None):
            self.phaseTransitionModel.setitem(self)
        if (self.relocationModel is not None):
            self.relocationModel.setitem(self)
        if (self.failureModel is not None):
            self.failureModel.setitem(self)
        if (self.poreVelocityModel is not None):
            self.poreVelocityModel.setitem(self)

        if (self.rheologyModel is not None):
            self.rheologyModel.setitem(self)

        return super().__repr__(depth)


    @property
    def Tref(self):
        return self._Tref

    @Tref.setter
    def Tref(self, Tref) -> None:
        check_type("Tref", Tref, (float, int))
        self._Tref = Tref

    @property
    def densityModel(self):
        return self._densityModel

    @densityModel.setter
    def densityModel(self, densityModel) -> None:
        check_type("densityModel", densityModel, DensityModel)
        self._densityModel = densityModel

    @property
    def heatCapacityModel(self):
        return self._heatCapacityModel

    @heatCapacityModel.setter
    def heatCapacityModel(self, heatCapacityModel) -> None:
        check_type("heatCapacityModel", heatCapacityModel, HeatCapacityModel)
        self._heatCapacityModel = heatCapacityModel

    @property
    def conductivityModel(self):
        return self._conductivityModel

    @conductivityModel.setter
    def conductivityModel(self, conductivityModel) -> None:
        check_type("conductivityModel", conductivityModel, ConductivityModel)
        self._conductivityModel = conductivityModel

    @property
    def emissivityModel(self):
        return self._emissivityModel

    @emissivityModel.setter
    def emissivityModel(self, emissivityModel) -> None:
        check_type("emissivityModel", emissivityModel, EmissivityModel)
        self._emissivityModel = emissivityModel

    @property
    def poissonRatioModel(self):
        return self._poissonRatioModel

    @poissonRatioModel.setter
    def poissonRatioModel(self, poissonRatioModel) -> None:
        check_type("poissonRatioModel", poissonRatioModel, PoissonRatioModel)
        self._poissonRatioModel = poissonRatioModel

    @property
    def thermalExpansionModel(self):
        return self._thermalExpansionModel

    @thermalExpansionModel.setter
    def thermalExpansionModel(self, thermalExpansionModel) -> None:
        check_type("thermalExpansionModel", thermalExpansionModel, ThermalExpansionModel)
        self._thermalExpansionModel = thermalExpansionModel

    @property
    def youngModulusModel(self):
        return self._youngModulusModel

    @youngModulusModel.setter
    def youngModulusModel(self, youngModulusModel) -> None:
        check_type("youngModulusModel", youngModulusModel, YoungModulusModel)
        self._youngModulusModel = youngModulusModel


    @property
    def densificationModel(self):
        return self._densificationModel

    @densificationModel.setter
    def densificationModel(self, densificationModel) -> None:
        check_type("densificationModel", densificationModel, DensificationModel, none_ok=True)
        self._densificationModel = densificationModel

    @property
    def swellingModel(self):
        return self._swellingModel

    @swellingModel.setter
    def swellingModel(self, swellingModel) -> None:
        check_type("swellingModel", swellingModel, SwellingModel, none_ok=True)
        self._swellingModel = swellingModel

    @property
    def phaseTransitionModel(self):
        return self._phaseTransitionModel

    @phaseTransitionModel.setter
    def phaseTransitionModel(self, phaseTransitionModel) -> None:
        check_type("phaseTransitionModel", phaseTransitionModel, PhaseTransitionModel, none_ok=True)
        self._phaseTransitionModel = phaseTransitionModel

    @property
    def relocationModel(self):
        return self._relocationModel

    @relocationModel.setter
    def relocationModel(self, relocationModel) -> None:
        check_type("relocationModel", relocationModel, RelocationModel, none_ok=True)
        self._relocationModel = relocationModel

    @property
    def failureModel(self):
        return self._failureModel

    @failureModel.setter
    def failureModel(self, failureModel) -> None:
        check_type("failureModel", failureModel, FailureModel, none_ok=True)
        self._failureModel = failureModel

    @property
    def poreVelocityModel(self):
        return self._poreVelocityModel

    @poreVelocityModel.setter
    def poreVelocityModel(self, poreVelocityModel) -> None:
        check_type("poreVelocityModel", poreVelocityModel, PoreVelocityModel, none_ok=True)
        self._poreVelocityModel = poreVelocityModel

    @property
    def rheologyModel(self):
        return self._rheologyModel

    @rheologyModel.setter
    def rheologyModel(self, rheologyModel) -> None:
        check_type("rheologyModel", rheologyModel, RheologyConstitutiveLaw, none_ok=True)
        self._rheologyModel = rheologyModel


class FuelMaterialModel(MaterialModel):
    def __init__(
            self,
            name: str,
            material: str,
            Tref: float,
            densityModel: DensityModel,
            heatCapacityModel: HeatCapacityModel,
            conductivityModel: ConductivityModel,
            emissivityModel: EmissivityModel,
            poissonRatioModel: PoissonRatioModel,
            thermalExpansionModel: ThermalExpansionModel,
            youngModulusModel: YoungModulusModel,
            densificationModel: DensificationModel,
            swellingModel: SwellingModel,
            relocationModel: RelocationModel,
            failureModel: FailureModel,
            poreVelocityModel: PoreVelocityModel,
            enrichment: float,
            rGrain: float,
            densityFraction: float,
            resinteringDensityChange: float,
            theoreticalDensity: float,
            nCracksMax: float,
            isotropicCracking: bool,
            GapCold: float,
            DiamCold: float,
            outerPatch: str,
            oxygenMetalRatio: float=2,
            ratioUMetal: float=0,
            ratioPuMetal: float=0,
            ratioAmMetal: float=0,
            dishFraction: float=0,
            burnupName: str="Bu",
            densityName: str="rho",
            heatSourceName: str="Q",
            gapWidthName: str="gapWidth",
            GdContent: float=0,
            nCracks: float=0,
            isotropicCrackingType: str="Barani",
            Tsintering: float=1800,
            recoveryFraction: float=0,
            modifiedRelocationModel: bool=True,
            rheologyModel: RheologyConstitutiveLaw=None
        ):
        super().__init__(
            name=name,
            material=material,
            Tref=Tref,
            densityModel=densityModel,
            heatCapacityModel=heatCapacityModel,
            conductivityModel=conductivityModel,
            emissivityModel=emissivityModel,
            poissonRatioModel=poissonRatioModel,
            thermalExpansionModel=thermalExpansionModel,
            youngModulusModel=youngModulusModel,
            densificationModel=densificationModel,
            swellingModel=swellingModel,
            relocationModel=relocationModel,
            failureModel=failureModel,
            poreVelocityModel=poreVelocityModel,
            rheologyModel=rheologyModel
        )

        self.enrichment = enrichment
        self.rGrain = rGrain
        self.densityFraction = densityFraction
        self.resinteringDensityChange = resinteringDensityChange
        self.theoreticalDensity = theoreticalDensity
        self.nCracksMax = nCracksMax
        self.isotropicCracking = isotropicCracking
        self.GapCold = GapCold
        self.DiamCold = DiamCold
        self.outerPatch = outerPatch
        self.oxygenMetalRatio = oxygenMetalRatio
        self.ratioUMetal = ratioUMetal
        self.ratioPuMetal = ratioPuMetal
        self.ratioAmMetal = ratioAmMetal
        self.dishFraction = dishFraction
        self.burnupName = burnupName
        self.densityName = densityName
        self.heatSourceName = heatSourceName
        self.gapWidthName = gapWidthName
        self.GdContent = GdContent
        self.nCracks = nCracks
        self.isotropicCrackingType = isotropicCrackingType
        self.Tsintering = Tsintering
        self.recoveryFraction = recoveryFraction
        self.modifiedRelocationModel = modifiedRelocationModel


    @property
    def enrichment(self):
        return self._enrichment

    @enrichment.setter
    def enrichment(self, enrichment) -> None:
        check_type("enrichment", enrichment, (float, int))
        check_positive("enrichment", enrichment)
        self._enrichment = enrichment
        self.__setitem__('enrichment', self.enrichment)

    @property
    def rGrain(self):
        return self._rGrain

    @rGrain.setter
    def rGrain(self, rGrain) -> None:
        check_type("rGrain", rGrain, (float, int))
        check_positive("rGrain", rGrain)
        self._rGrain = rGrain
        self.__setitem__('rGrain', self.rGrain)

    @property
    def densityFraction(self):
        return self._densityFraction

    @densityFraction.setter
    def densityFraction(self, densityFraction) -> None:
        check_type("densityFraction", densityFraction, (float, int))
        check_positive("densityFraction", densityFraction)
        self._densityFraction = densityFraction
        self.__setitem__('densityFraction', self.densityFraction)

    @property
    def resinteringDensityChange(self):
        return self._resinteringDensityChange

    @resinteringDensityChange.setter
    def resinteringDensityChange(self, resinteringDensityChange) -> None:
        check_type("resinteringDensityChange", resinteringDensityChange, (float, int))
        self._resinteringDensityChange = resinteringDensityChange
        self.__setitem__('resinteringDensityChange', self.resinteringDensityChange)

    @property
    def theoreticalDensity(self):
        return self._theoreticalDensity

    @theoreticalDensity.setter
    def theoreticalDensity(self, theoreticalDensity) -> None:
        check_type("theoreticalDensity", theoreticalDensity, (float, int))
        check_positive("theoreticalDensity", theoreticalDensity)
        self._theoreticalDensity = theoreticalDensity
        self.__setitem__('theoreticalDensity', self.theoreticalDensity)

    @property
    def nCracksMax(self):
        return self._nCracksMax

    @nCracksMax.setter
    def nCracksMax(self, nCracksMax) -> None:
        check_type("nCracksMax", nCracksMax, (float, int))
        check_positive("nCracksMax", nCracksMax)
        self._nCracksMax = nCracksMax
        self.__setitem__('nCracksMax', self.nCracksMax)

    @property
    def isotropicCracking(self):
        return self._isotropicCracking

    @isotropicCracking.setter
    def isotropicCracking(self, isotropicCracking) -> None:
        check_type("isotropicCracking", isotropicCracking, bool)
        self._isotropicCracking = isotropicCracking
        self.__setitem__('isotropicCracking', self.isotropicCracking)

    @property
    def GapCold(self):
        return self._GapCold

    @GapCold.setter
    def GapCold(self, GapCold) -> None:
        check_type("GapCold", GapCold, (float, int))
        self._GapCold = GapCold
        self.__setitem__('GapCold', self.GapCold)

    @property
    def DiamCold(self):
        return self._DiamCold

    @DiamCold.setter
    def DiamCold(self, DiamCold) -> None:
        check_type("DiamCold", DiamCold, (float, int))
        check_positive("DiamCold", DiamCold)
        self._DiamCold = DiamCold
        self.__setitem__('DiamCold', self.DiamCold)

    @property
    def outerPatch(self):
        return self._outerPatch

    @outerPatch.setter
    def outerPatch(self, outerPatch) -> None:
        check_type("outerPatch", outerPatch, str)
        self._outerPatch = outerPatch
        self.__setitem__('outerPatch', self.outerPatch)

    @property
    def oxygenMetalRatio(self):
        return self._oxygenMetalRatio

    @oxygenMetalRatio.setter
    def oxygenMetalRatio(self, oxygenMetalRatio) -> None:
        check_type("oxygenMetalRatio", oxygenMetalRatio, (float, int))
        check_positive("oxygenMetalRatio", oxygenMetalRatio)
        self._oxygenMetalRatio = oxygenMetalRatio
        self.__setitem__('oxygenMetalRatio', self.oxygenMetalRatio)

    @property
    def ratioUMetal(self):
        return self._ratioUMetal

    @ratioUMetal.setter
    def ratioUMetal(self, ratioUMetal) -> None:
        check_type("ratioUMetal", ratioUMetal, (float, int))
        check_positive("ratioUMetal", ratioUMetal)
        self._ratioUMetal = ratioUMetal
        self.__setitem__('ratioUMetal', self.ratioUMetal)

    @property
    def ratioPuMetal(self):
        return self._ratioPuMetal

    @ratioPuMetal.setter
    def ratioPuMetal(self, ratioPuMetal) -> None:
        check_type("ratioPuMetal", ratioPuMetal, (float, int))
        check_positive("ratioPuMetal", ratioPuMetal)
        self._ratioPuMetal = ratioPuMetal
        self.__setitem__('ratioPuMetal', self.ratioPuMetal)

    @property
    def ratioAmMetal(self):
        return self._ratioAmMetal

    @ratioAmMetal.setter
    def ratioAmMetal(self, ratioAmMetal) -> None:
        check_type("ratioAmMetal", ratioAmMetal, (float, int))
        check_positive("ratioAmMetal", ratioAmMetal)
        self._ratioAmMetal = ratioAmMetal
        self.__setitem__('ratioAmMetal', self.ratioAmMetal)

    @property
    def dishFraction(self):
        return self._dishFraction

    @dishFraction.setter
    def dishFraction(self, dishFraction) -> None:
        check_type("dishFraction", dishFraction, (float, int))
        self._dishFraction = dishFraction
        self.__setitem__('dishFraction', self.dishFraction)

    @property
    def burnupName(self):
        return self._burnupName

    @burnupName.setter
    def burnupName(self, burnupName) -> None:
        check_type("burnupName", burnupName, str)
        self._burnupName = burnupName
        self.__setitem__('burnupName', self.burnupName)

    @property
    def densityName(self):
        return self._densityName

    @densityName.setter
    def densityName(self, densityName) -> None:
        check_type("densityName", densityName, str)
        self._densityName = densityName
        self.__setitem__('densityName', self.densityName)

    @property
    def heatSourceName(self):
        return self._heatSourceName

    @heatSourceName.setter
    def heatSourceName(self, heatSourceName) -> None:
        check_type("heatSourceName", heatSourceName, str)
        self._heatSourceName = heatSourceName
        self.__setitem__('heatSourceName', self.heatSourceName)

    @property
    def gapWidthName(self):
        return self._gapWidthName

    @gapWidthName.setter
    def gapWidthName(self, gapWidthName) -> None:
        check_type("gapWidthName", gapWidthName, str)
        self._gapWidthName = gapWidthName
        self.__setitem__('gapWidthName', self.gapWidthName)

    @property
    def GdContent(self):
        return self._GdContent

    @GdContent.setter
    def GdContent(self, GdContent) -> None:
        check_type("GdContent", GdContent, (float, int))
        check_positive("GdContent", GdContent)
        self._GdContent = GdContent
        self.__setitem__('GdContent', self.GdContent)

    @property
    def nCracks(self):
        return self._nCracks

    @nCracks.setter
    def nCracks(self, nCracks) -> None:
        check_type("nCracks", nCracks, (float, int))
        check_positive("nCracks", nCracks)
        self._nCracks = nCracks
        self.__setitem__('nCracks', self.nCracks)

    @property
    def isotropicCrackingType(self):
        return self._isotropicCrackingType

    @isotropicCrackingType.setter
    def isotropicCrackingType(self, isotropicCrackingType) -> None:
        check_type("isotropicCrackingType", isotropicCrackingType, str)
        check_value("isotropicCrackingType", isotropicCrackingType, {"Barani", "JankusWeeks"})
        self._isotropicCrackingType = isotropicCrackingType
        self.__setitem__('isotropicCrackingType', self.isotropicCrackingType)

    @property
    def Tsintering(self):
        return self._Tsintering

    @Tsintering.setter
    def Tsintering(self, Tsintering) -> None:
        check_type("Tsintering", Tsintering, (float, int))
        check_positive("Tsintering", Tsintering)
        self._Tsintering = Tsintering
        self.__setitem__('Tsintering', self.Tsintering)

    @property
    def recoveryFraction(self):
        return self._recoveryFraction

    @recoveryFraction.setter
    def recoveryFraction(self, recoveryFraction) -> None:
        check_type("recoveryFraction", recoveryFraction, (float, int))
        check_positive("recoveryFraction", recoveryFraction)
        self._recoveryFraction = recoveryFraction
        self.__setitem__('recoveryFraction', self.recoveryFraction)

    @property
    def modifiedRelocationModel(self):
        return self._modifiedRelocationModel

    @modifiedRelocationModel.setter
    def modifiedRelocationModel(self, modifiedRelocationModel) -> None:
        check_type("modifiedRelocationModel", modifiedRelocationModel, bool)
        self._modifiedRelocationModel = modifiedRelocationModel
        self.__setitem__('modifiedRelocationModel', self.modifiedRelocationModel)


class ConstantMechanicalLaw:
    def __init__(
            self,
            name: str,
            dimension: ReducedDimension,
            value: float
        ):
        self.name = name
        self.dimension = dimension
        self.value = value

    def __repr__(self):
        return(f"{self.name:11} {self.dimension!r:15} {self.value:g}")


    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name) -> None:
        check_type("name", name, str)
        self._name = name

    @property
    def dimension(self):
        return self._dimension

    @dimension.setter
    def dimension(self, dimension) -> None:
        check_type("dimension", dimension, ReducedDimension)
        self._dimension = dimension

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value) -> None:
        check_type("value", value, (float, int))
        self._value = value
