from foamForNuclear.common import FoamForNuclearDict
from foamForNuclear._attrs_tools import ffn_define, _to_List_float
from typing import ClassVar
from foamForNuclear.common import Table, List
from attrs import field


@ffn_define
class Creep(FoamForNuclearDict):
    TYPE: ClassVar[str] = "none"


@ffn_define
class ConstantPrincipalStress(Creep):
    TYPE: ClassVar[str] = "constantPrincipalStress"
    relax: float | int = 1.0
    fastFluxName: str = "fastFlux"
    poissonRationForCreep: float | int = 0.5
    fluxConversionFactor: float | int = 1.0
    creepCoefficient: float | int


@ffn_define
class CorrelationPrincipalStress(Creep):
    TYPE: ClassVar[str] = "correlationPrincipalStress"
    relax: float | int = 1.0
    fastFluxName: str = "fastFlux"
    poissonRationForCreep: float | int = 0.5
    fluxConversionFactor: float | int = 1.0
    coefficientList: list | List = field(factory=list, converter=_to_List_float)


@ffn_define
class Limback(Creep):
    TYPE: ClassVar[str] = "ZircaloyLimback"
    relax: float | int = 1.0
    fastFluxName: str = "fastFlux"
    fastFluenceName: str = "fastFluence"
    cladType: str = "SRA"


@ffn_define
class LimbackLoca(Creep):
    TYPE: ClassVar[str] = "ZircaloyLimbackLoca"
    relax: float | int = 1.0
    fastFluxName: str = "fastFlux"
    fastFluenceName: str = "fastFluence"
    irradiationCreep: bool = True
    primaryCreep: bool = True
    primaryCreepAbove700K: bool = True
    LOCACreep: bool = True
    NewtonRaphsonMethod: bool = False
    NewtonRaphsonTolerance: float | int = 1e-6


@ffn_define
class MOXMalygin(Creep):
    TYPE: ClassVar[str] = "UPuO2Malygin"
    relax: float | int = 1.0
    heatSourceName: str = "Q"
    grainRadiusName: str = "grainRadius"
    oxygenMetalRatioName: str = "oxygenMetalRatio"
    youngModulusName: str = "E"
    porosityName: str = "porosity"
    PuConcentration: float | int
    irradiationCreep: bool = True
    NewtonRaphsonMethod: bool = False
    NewtonRaphsonTolerance: float | int = 1e-6
    creepRateLimiter: bool = False
    creepRateLimit: float | int = 4.1e-6


@ffn_define
class UO2Matpro(Creep):
    TYPE: ClassVar[str] = "UO2Matpro"
    relax: float | int = 1.0
    heatSourceName: str = "Q"
    densityName: str = "rho"
    grainRadiusName: str = "grainRadius"
    oxygenMetalRatioName: str = "oxygenMetalRatio"
    SakaiCorrection: bool = True
    NewtonRaphsonMethod: bool = False
    NewtonRaphsonTolerance: float | int = 1e-6


@ffn_define
class SiCMonolithic(Creep):
    TYPE: ClassVar[str] = "SiCMonolithic"
    relax: float | int = 1.0
    fastFluxName: str = "fastFlux"
    creepCoefficient: float | int = -1.0
    poissonRationForCreep: float | int = 0.5
    fluxConversionFactor: float | int = 1.0


@ffn_define
class BufferParfume(Creep):
    TYPE: ClassVar[str] = "BufferParfume"
    relax: float | int = 1.0
    densityName: str = "rho"
    fastFluxName: str = "fastFlux"
    creepAmplificationCoefficient: float | int = 2.0
    poissonRationForCreep: float | int = 0.5
    fluxConversionFactor: float | int = 1.0


@ffn_define
class PyCParfume(Creep):
    TYPE: ClassVar[str] = "PyCParfume"
    relax: float | int = 1.0
    densityName: str = "rho"
    fastFluxName: str = "fastFlux"
    creepAmplificationCoefficient: float | int = 2.0
    poissonRationForCreep: float | int = 0.5
    fluxConversionFactor: float | int = 1.0


@ffn_define
class PowerLaw(Creep):
    TYPE: ClassVar[str] = "powerLaw"
    relax: float | int = 1.0
    B: float | int
    sigmaC: float | int
    n: float | int


@ffn_define
class MOXRoutbort(Creep):
    TYPE: ClassVar[str] = "UPuO2Routbort"
    relax: float | int = 1.0
    heatSourceName: str = "Q"
    grainRadiusName: str = "grainRadius"
    irradiationCreep: bool = True
    NewtonRaphsonMethod: bool = False
    NewtonRaphsonTolerance: float | int = 1e-6
    creepRateLimiter: bool = False
    creepRateLimit: float | int = 4.1e-6


@ffn_define
class Steel1515TiAIM1Tobbe(Creep):
    TYPE: ClassVar[str] = "Steel1515TiAim1Tobbe"
    relax: float | int = 1.0
    fastFluxName: str = "fastFlux"


@ffn_define
class Steel1515TiTobbe(Creep):
    TYPE: ClassVar[str] = "Steel1515TiTobbe"
    relax: float | int = 1.0
    fastFluxName: str = "fastFlux"
    fastFluenceName: str = "fastFluence"


@ffn_define
class Steel1515TiDINTobbe(Creep):
    TYPE: ClassVar[str] = "Steel1515TiDinTobbe"
    relax: float | int = 1.0
    fastFluxName: str = "fastFlux"


@ffn_define
class HastelloyZhang(Creep):
    TYPE: ClassVar[str] = "HastelloyZhang"
    relax: float | int = 1.0
    fastFluxName: str = "fastFlux"
    fastFluenceName: str = "fastFluence"    
    NewtonRaphsonMethod: bool = False
    NewtonRaphsonTolerance: float | int = 1e-6