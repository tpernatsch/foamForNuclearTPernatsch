from foamForNuclear.thermomechanicalMaterial.thermomechanicalMaterial import *


class DensityIAEAZy(DensityModel):
    def __init__(
            self,
            par1: float=6595.2,
            par2: float=0.1477,
            par3: float=6690.0,
            par4: float=0.1855,
            perturb: float=1.0
        ):
        super().__init__(name="ZyIAEA")

        self.par1 = par1
        self.par2 = par2
        self.par3 = par3
        self.par4 = par4
        self.perturb = perturb


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("density", OpenFOAMDict({
            "par1": self.par1,
            "par2": self.par2,
            "par3": self.par3,
            "par4": self.par4,
            "perturb": self.perturb,
        }))


class HeatCapacityIAEAZy(HeatCapacityModel):
    def __init__(
            self,
            par1: float=255.66,
            par2: float=0.1024,
            par3: float=597.1,
            par4: float=0.4088,
            par5: float=1.565E-4,
            par6: float=1058.4,
            par7: float=1213.8,
            par8: float=2.0,
            par9: float=719.61,
            perturb: float=1.0
        ):
        super().__init__(name="ZyIAEA")

        self.par1 = par1
        self.par2 = par2
        self.par3 = par3
        self.par4 = par4
        self.par5 = par5
        self.par6 = par6
        self.par7 = par7
        self.par8 = par8
        self.par9 = par9
        self.perturb = perturb


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("heatCapacity", OpenFOAMDict({
            "par1": self.par1,
            "par2": self.par2,
            "par3": self.par3,
            "par4": self.par4,
            "par5": self.par5,
            "par6": self.par6,
            "par7": self.par7,
            "par8": self.par8,
            "par9": self.par9,
            "perturb": self.perturb,
        }))


class ConductivityRelapZy(ConductivityModel):
    def __init__(
            self,
            par1: float=7.51,
            par2: float=2.09e-2,
            par3: float=1.45e-5,
            par4: float=7.67e-9,
            perturb: float=1.0
        ):
        super().__init__(name="ZyRELAP")

        self.par1 = par1
        self.par2 = par2
        self.par3 = par3
        self.par4 = par4
        self.perturb = perturb


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("conductivity", OpenFOAMDict({
            "par1": self.par1,
            "par2": self.par2,
            "par3": self.par3,
            "par4": self.par4,
            "perturb": self.perturb,
        }))


class ConstantEmissivityZy(EmissivityModel):
    def __init__(
            self,
            emissivityValue: float=0.808642
        ):
        super().__init__(name="ZyConstant")

        self.emissivityValue = emissivityValue


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("emissivity", OpenFOAMDict({
            "emissivityValue": self.emissivityValue,
        }))


class YoungModulusMatproZy(YoungModulusModel):
    def __init__(
            self,
            coldWork: float=0,
            oxyCon: float=0,
            par1: float=6.61e11,
            par2: float=5.912e8,
            par3: float=2.6e10,
            par4: float=0.88,
            par5: float=0.12,
            par6: float=1e25,
            par7: float=1.088e11,
            par8: float=5.475e7,
            par9: float=9.21e10,
            par10: float=4.05e7,
            perturb: float=1.0
        ):
        super().__init__(name="ZyMATPRO")

        self.coldWork = coldWork
        self.oxyCon = oxyCon

        self.par1 = par1
        self.par2 = par2
        self.par3 = par3
        self.par4 = par4
        self.par5 = par5
        self.par6 = par6
        self.par7 = par7
        self.par8 = par8
        self.par9 = par9
        self.par10 = par10
        self.perturb = perturb


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("coldWork", self.coldWork)
        parentMaterial.__setitem__("oxyCon", self.oxyCon)
        parentMaterial.__setitem__("YoungModulus", OpenFOAMDict({
            "par1": self.par1,
            "par2": self.par2,
            "par3": self.par3,
            "par4": self.par4,
            "par5": self.par5,
            "par6": self.par6,
            "par7": self.par7,
            "par8": self.par8,
            "par9": self.par9,
            "par10": self.par10,
            "perturb": self.perturb,
        }))


class ConstantPoissonRatioZy(PoissonRatioModel):
    def __init__(
            self,
            poissonRatioValue: float=0.3
        ):
        super().__init__(name="ZyConstant")

        self.poissonRatioValue = poissonRatioValue


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("PoissonRatio", OpenFOAMDict({
            "PoissonRatioValue": self.poissonRatioValue,
        }))


class ThermalExpansionMatproZy(ThermalExpansionModel):
    def __init__(
            self,
            par1: float=4.441e-6,
            par2: float=1.238e-3,
            par3: float=6.721e-6,
            par4: float=2.073e-3,
            par5: float=9.7e-6,
            par6: float=1.10e-2,
            par7: float=9.45e-3,
            perturb: float=1.0
        ):
        super().__init__(name="ZyMATPRO")

        self.par1 = par1
        self.par2 = par2
        self.par3 = par3
        self.par4 = par4
        self.par5 = par5
        self.par6 = par6
        self.par7 = par7
        self.perturb = perturb


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("thermalExpansion", OpenFOAMDict({
            "par1": self.par1,
            "par2": self.par2,
            "par3": self.par3,
            "par4": self.par4,
            "par5": self.par5,
            "par6": self.par6,
            "par7": self.par7,
            "perturb": self.perturb,
        }))


class SwellingGrowthBISONZy(SwellingModel):
    def __init__(
            self,
            A: float=3e-20,
            n: float=0.794,
            perturb: float=1.0
        ):
        super().__init__(name="growthBISONZy")

        self.A = A
        self.n = n
        self.perturb = perturb


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("swelling", OpenFOAMDict({
            "A": self.A,
            "n": self.n,
            "perturb": self.perturb,
        }))


class PhaseTransitionZyDynamic(PhaseTransitionModel):
    def __init__(self):
        super().__init__(name="ZyDynamic")


class Zircalloy(MaterialModel):
    def __init__(
        self,
        name,
        Tref,
        densityModel: DensityModel=None,
        heatCapacityModel: HeatCapacityModel=None,
        conductivityModel: ConductivityModel=None,
        emissivityModel: EmissivityModel=None,
        poissonRatioModel: PoissonRatioModel=None,
        thermalExpansionModel: ThermalExpansionModel=None,
        youngModulusModel: YoungModulusModel=None,
        swellingModel: SwellingModel=None,
        phaseTransitionModel: PhaseTransitionModel=None,
        failureModel: FailureModel=None,
        rheologyModel: RheologyConstitutiveLaw=None
    ):
        super().__init__(
            name=name,
            material="zircaloy",
            Tref=Tref,
            densityModel=DensityIAEAZy() if densityModel is None else densityModel,
            heatCapacityModel=HeatCapacityIAEAZy() if heatCapacityModel is None else heatCapacityModel,
            conductivityModel=ConductivityRelapZy() if conductivityModel is None else conductivityModel,
            emissivityModel=ConstantEmissivityZy() if emissivityModel is None else emissivityModel,
            poissonRatioModel=ConstantPoissonRatioZy() if poissonRatioModel is None else poissonRatioModel,
            thermalExpansionModel=ThermalExpansionMatproZy() if thermalExpansionModel is None else thermalExpansionModel,
            youngModulusModel=YoungModulusMatproZy() if youngModulusModel is None else youngModulusModel,
            swellingModel=SwellingGrowthBISONZy() if swellingModel is None else swellingModel,
            phaseTransitionModel=PhaseTransitionZyDynamic() if phaseTransitionModel is None else phaseTransitionModel,
            failureModel=FailureModel() if failureModel is None else failureModel,
            rheologyModel=rheologyModel
        )
