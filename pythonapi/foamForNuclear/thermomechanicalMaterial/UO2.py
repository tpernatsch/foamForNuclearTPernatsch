

import numpy as np
from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import OpenFOAMDict, addParameter
from foamForNuclear.thermomechanicalMaterial.thermomechanicalMaterial import ConductivityModel, DensificationModel, DensityModel, EmissivityModel, FailureModel, FuelMaterialModel, HeatCapacityModel, MaterialModel, PoissonRatioModel, PoreVelocityModel, RelocationModel, RheologyConstitutiveLaw, SwellingModel, ThermalExpansionModel, YoungModulusModel


class ConductivityMatproUO2(ConductivityModel):
    def __init__(
            self,
            par1: float=0.0452,
            par2: float=0.000246,
            par3: float=0.00187,
            par4: float=1.1599,
            par5: float=0.9,
            par6: float=0.04,
            par7: float=0.038,
            par8: float=0.28,
            par9: float=396.0,
            par10: float=6380.0,
            par11: float=3.5e9,
            par12: float=2.0,
            par13: float=16360.0,
            perturb: float=1.0
        ):
        super().__init__(name="UO2MATPRO")

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
        self.par11 = par11
        self.par12 = par12
        self.par13 = par13
        self.perturb = perturb


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("conductivity", OpenFOAMDict({
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
            "par11": self.par11,
            "par12": self.par12,
            "par13": self.par13,
            "perturb": self.perturb,
        }))


class ConstantDensityUO2(DensityModel):
    """
    Class to set UO2 density to a constant input value.
    """
    def __init__(self):
        super().__init__(name="UO2Constant")


class EmissivityRelapUO2(EmissivityModel):
    def __init__(
            self,
            par1: float=0.7856,
            par2: float=1.5263e-5,
            perturb: float=1.0
        ):
        super().__init__(name="UO2RELAP")

        self.par1 = par1
        self.par2 = par2
        self.perturb = perturb


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("emissivity", OpenFOAMDict({
            "par1": self.par1,
            "par2": self.par2,
            "perturb": self.perturb,
        }))


class HeatCapacityMatproUO2(HeatCapacityModel):
    """
    Heat Capacity model for UO2 derived from MATPROv11.
    """
    def __init__(
            self,
            OM: float=2.0,
            K1: float=296.7,
            K2: float=2.43e-2,
            K3: float=8.745e7,
            phi: float=535.285,
            ED: float=1.577e5,
            par1: float=2.0,
            par2: float=2.0,
            par3: float=2.0,
            perturb: float=1.0,
        ):
        super().__init__(name="UO2MATPRO")

        self.OM = OM
        self.K1 = K1
        self.K2 = K2
        self.K3 = K3
        self.phi = phi
        self.ED = ED
        self.par1 = par1
        self.par2 = par2
        self.par3 = par3
        self.perturb = perturb


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("heatCapacity", OpenFOAMDict({
            "OM": self.OM,
            "K1": self.K1,
            "K2": self.K2,
            "K3": self.K3,
            "phi": self.phi,
            "ED": self.ED,
            "par1": self.par1,
            "par2": self.par2,
            "par3": self.par3,
            "perturb": self.perturb,
        }))

    def value(self, T: float):
        R = 1.380649e-23 * 6.02214076e23

        nominalValue = (
            self.K1 * np.pow(self.phi, self.par1) * np.exp(self.phi/T) / np.pow(T*(np.exp(self.phi/T)-1), self.par2)
            + self.K2 * T
            + (self.OM/2) * self.K3 * self.ED/(R*np.pow(T, self.par3))*np.exp(-self.ED/(R*T))
        )

        return(self.perturb * nominalValue)


class ConstantPoissonRatioUO2(PoissonRatioModel):
    def __init__(
            self,
            PoissonRatioValue: float=0.316
        ):
        super().__init__(name="UO2Constant")

        self.PoissonRatioValue = PoissonRatioValue


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("PoissonRatio", OpenFOAMDict({
            "PoissonRatioValue": self.PoissonRatioValue,
        }))


class ThermalExpansionRelapUO2(ThermalExpansionModel):
    def __init__(
            self,
            K1: float=9.8e-6,
            K2: float=2.61e-3,
            K3: float=3.16e-1,
            ED: float=1.32e-19,
            k: float=1.38e-23,
            perturb: float=1.0
        ):
        super().__init__(name="UO2RELAP")

        self.K1 = K1
        self.K2 = K2
        self.K3 = K3
        self.ED = ED
        self.k = k
        self.perturb = perturb


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("thermalExpansion", OpenFOAMDict({
            "K1": self.K1,
            "K2": self.K2,
            "K3": self.K3,
            "ED": self.ED,
            "k": self.k,
            "perturb": self.perturb,
        }))


class YoungModulusMatproUO2(YoungModulusModel):
    def __init__(
            self,
            par1: float=2.334e11,
            par2: float=2.752,
            par3: float=1.0915e-4,
            perturb: float=1.0
        ):
        super().__init__(name="UO2MATPRO")

        self.par1 = par1
        self.par2 = par2
        self.par3 = par3
        self.perturb = perturb


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("YoungModulus", OpenFOAMDict({
            "par1": self.par1,
            "par2": self.par2,
            "par3": self.par3,
            "perturb": self.perturb,
        }))


class DensificationFrapcon(DensificationModel):
    def __init__(
            self,
            par1: float=22.2,
            par2: float=1453,
            par3: float=66.6,
            par4: float=2.0,
            par5: float=35.0,
            par6: float=3.0,
            perturb: float=1.0,
        ):
        super().__init__(name="UO2FRAPCON")

        self.par1 = par1
        self.par2 = par2
        self.par3 = par3
        self.par4 = par4
        self.par5 = par5
        self.par6 = par6
        self.perturb = perturb


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("densification", OpenFOAMDict({
            "par1": self.par1,
            "par2": self.par2,
            "par3": self.par3,
            "par4": self.par4,
            "par5": self.par5,
            "par6": self.par6,
            "perturb": self.perturb,
        }))


class SwellingFrapcon(SwellingModel):
    def __init__(
            self,
            par1: float=6000,
            par2: float=80000,
            par3: float=2.974e+10,
            par4: float=2.315e-23,
            par5: float=86.4,
            par6: float=3.211e-23,
            perturb: float=1.0,
        ):
        super().__init__(name="UO2FRAPCON")

        self.par1 = par1
        self.par2 = par2
        self.par3 = par3
        self.par4 = par4
        self.par5 = par5
        self.par6 = par6
        self.perturb = perturb


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("swelling", OpenFOAMDict({
            "par1": self.par1,
            "par2": self.par2,
            "par3": self.par3,
            "par4": self.par4,
            "par5": self.par5,
            "par6": self.par6,
            "perturb": self.perturb,
        }))


class RelocationFrapcon(RelocationModel):
    def __init__(
            self,
            perturb: float=1.0
        ):
        super().__init__(name="UO2FRAPCON")

        self.perturb = perturb


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("relocation", OpenFOAMDict({
            "perturb": self.perturb,
        }))


class FailureUO2meltingMatpro(FailureModel):
    def __init__(self, Tmelt: float=3113):
        super().__init__(name="UO2meltingMATPRO")

        self.Tmelt = Tmelt


    def setitem(self, parentMaterial: OpenFOAMDict):
        super().setitem(parentMaterial)
        parentMaterial.__setitem__("Tmelt", self.Tmelt)


class PoreVelocityUO2Sens(PoreVelocityModel):
    def __init__(self):
        super().__init__(name="UO2Sens")


class UO2(FuelMaterialModel):
    """
    UO2 material object.

    Parameters
    ----------
    densityModel: DensityModel
        Default `ConstantDensityUO2()`
    heatCapacityModel: HeatCapacityModel
        Default `HeatCapacityMatproUO2()`
    conductivityModel: ConductivityModel
        Default `ConductivityMatproUO2()`
    emissivityModel : EmissivityModel
        Default `EmissivityRelapUO2()`
    poissonRatioModel : PoissonRatioModel
        Default `ConstantPoissonRatioUO2()`
    thermalExpansionModel : ThermalExpansionModel
        Default `ThermalExpansionRelapUO2()`
    youngModulusModel : YoungModulusModel
        Default `YoungModulusMatproUO2()`
    """
    def __init__(
        self,
        name: str,
        Tref: float,
        enrichment: float,
        rGrain: float,
        nCracksMax: float,
        isotropicCracking: bool,
        GapCold: float,
        DiamCold: float,
        outerPatch: str,
        resinteringDensityChange: float,
        densityFraction: float=0.945,
        theoreticalDensity: float=10430.0,
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
        densityModel: DensityModel=None,
        heatCapacityModel: HeatCapacityModel=None,
        conductivityModel: ConductivityModel=None,
        emissivityModel: EmissivityModel=None,
        poissonRatioModel: PoissonRatioModel=None,
        thermalExpansionModel: ThermalExpansionModel=None,
        youngModulusModel: YoungModulusModel=None,
        densificationModel: DensificationModel=None,
        swellingModel: SwellingModel=None,
        relocationModel: RelocationModel=None,
        failureModel: FailureModel=None,
        poreVelocityModel: PoreVelocityModel=None,
        rheologyModel: RheologyConstitutiveLaw=None
    ):
        super().__init__(
            name=name,
            material="UO2",
            Tref=Tref,
            densityModel=ConstantDensityUO2() if densityModel is None else densityModel,
            heatCapacityModel=HeatCapacityMatproUO2() if heatCapacityModel is None else heatCapacityModel,
            conductivityModel=ConductivityMatproUO2() if conductivityModel is None else conductivityModel,
            emissivityModel=EmissivityRelapUO2() if emissivityModel is None else emissivityModel,
            poissonRatioModel=ConstantPoissonRatioUO2() if poissonRatioModel is None else poissonRatioModel,
            thermalExpansionModel=ThermalExpansionRelapUO2() if thermalExpansionModel is None else thermalExpansionModel,
            youngModulusModel=YoungModulusMatproUO2() if youngModulusModel is None else youngModulusModel,
            densificationModel=DensificationFrapcon() if densificationModel is None else densificationModel,
            swellingModel=SwellingFrapcon() if swellingModel is None else swellingModel,
            relocationModel=RelocationFrapcon() if relocationModel is None else relocationModel,
            failureModel=FailureModel() if failureModel is None else failureModel,
            poreVelocityModel=PoreVelocityUO2Sens() if poreVelocityModel is None else poreVelocityModel,
            rheologyModel=rheologyModel,
            enrichment=enrichment,
            rGrain=rGrain,
            densityFraction=densityFraction,
            resinteringDensityChange=resinteringDensityChange,
            theoreticalDensity=theoreticalDensity,
            nCracksMax=nCracksMax,
            isotropicCracking=isotropicCracking,
            GapCold=GapCold,
            DiamCold=DiamCold,
            outerPatch=outerPatch,
            oxygenMetalRatio=oxygenMetalRatio,
            ratioUMetal=ratioUMetal,
            ratioPuMetal=ratioPuMetal,
            ratioAmMetal=ratioAmMetal,
            dishFraction=dishFraction,
            burnupName=burnupName,
            densityName=densityName,
            heatSourceName=heatSourceName,
            gapWidthName=gapWidthName,
            GdContent=GdContent,
            nCracks=nCracks,
            isotropicCrackingType=isotropicCrackingType,
            Tsintering=Tsintering,
            recoveryFraction=recoveryFraction,
            modifiedRelocationModel=modifiedRelocationModel,
        )
