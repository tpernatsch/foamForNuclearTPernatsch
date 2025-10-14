from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import OpenFOAMDict


_HEAT_TRANSFER_MODEL_TYPES = {
    "CachardLiquid", "CachardVapour", "Gorenflo", "NusseltAndWall",
    "NusseltReynoldsPrandtlPower", "NusseltWallAndHfromFMU", "Shah", "byRegime",
    "constant", "multiRegimeBoiling", "multiRegimeBoilingTRACE",
    "multiRegimeBoilingTRACECHF", "multiRegimeBoilingVapourTRACE",
    "superpositionNucleateBoiling",
    "NoKazimi"
}

class HeatTransferModel(OpenFOAMDict):
    """
    Base class for heat-transfer.
    """
    def __init__(
            self,
            type: str,
            zones: list[str]=[],
        ):
        super().__init__()
        self.type = type
        self.zones = zones


    def __repr__(self, depth: int=0):
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
        check_value("type", type, _HEAT_TRANSFER_MODEL_TYPES)
        self._type_ = type
        self.__setitem__("type", type)


class Gorenflo(HeatTransferModel):
    r"""
    Gorenflo Heat Transfer Model

    Implementation Notes
    --------------------
    The original pool boiling heat transfer coefficient by Gorenflo is the following
    form:

    .. math::
        htcPB = h_0 F ((\frac{q}{q_0})^n)*(\frac{R}{R_0})^0.133

    with F being a function of the reduced pressure, i.e. p/pCtri. Needless to
    say, this form is not particularly suitable for a numerical implementation in
    its current form as htcPB is required to compute q in the first place. In
    general, the pool boiling heat flux is given by:

    .. math::
        q = htc*(Twall-Tsat)

    with htc being the total heat transfer coefficient between wall and fluid,
    Twall being the wall temperature and Tsat the saturation temperature. By
    assuming that htc ~ htcPB (i.e. that most of the heat flux is due to the pool
    boiling mechanism), one can tha equate:

    .. math::
        htcPB = h_0 F ((htcPB*(Twall-Tsat)/q_0)^n)*(\frac{R}{R_0})^0.133

    which can be re-arranged to isolate htcPB on the left hand side as:

    .. math::
        htcPB = (h_0 F ((Twall-Tsat)^n)*(\frac{R}{R_0})^0.133/(q_0^n))^(1.0/(1.0-n))

    Which is how the correlation is implemented here.
    The assumption of htc ~ htcPB can be rather crude but it is widely employed
    to yield a feasable numerical implementation in other computer codes too (e.g.
    TRACE)


    Parameters
    ----------
    absoluteSurfaceRoughness : float
        Default `4e-7`
    useExplicitHeatFlux : bool
        Default `False`
    """

    def __init__(
            self,
            absoluteSurfaceRoughness: float=4e-7,
            useExplicitHeatFlux: bool=False,
            zones=[],
        ):
        super().__init__("Gorenflo", zones)

        self.absoluteSurfaceRoughness = absoluteSurfaceRoughness
        self.useExplicitHeatFlux = useExplicitHeatFlux


    def __repr__(self, depth=0):
        self.__setitem__("absoluteSurfaceRoughness", self.absoluteSurfaceRoughness)
        self.__setitem__("useExplicitHeatFlux", self.useExplicitHeatFlux)

        return super().__repr__(depth)


class NusseltReynoldsPrandtlPower(HeatTransferModel):
    """
    Nusselt Heat Transfer Model

    Parameters
    ----------
    const : float
    coeff : float
    expRe : float
    expPr : float
    expTc : float
        Default 0
    """

    def __init__(
            self,
            const: float,
            coeff: float,
            expRe: float,
            expPr: float,
            expTc: float=0,
            zones=[],
        ):
        super().__init__("NusseltReynoldsPrandtlPower", zones)

        self.const = const
        self.coeff = coeff
        self.expRe = expRe
        self.expPr = expPr
        self.expTc = expTc


    def __repr__(self, depth=0):
        self.__setitem__("const", self.const)
        self.__setitem__("coeff", self.coeff)
        self.__setitem__("expRe", self.expRe)
        self.__setitem__("expPr", self.expPr)
        self.__setitem__("expTc", self.expTc)

        return super().__repr__(depth)


class NusseltAndWall(HeatTransferModel):
    """
    Nusselt Heat Transfer Model with an additional wall

    Parameters
    ----------
    const : float
    coeff : float
    expRe : float
    expPr : float
    addH : float
    """

    def __init__(
            self,
            const: float,
            coeff: float,
            expRe: float,
            expPr: float,
            addH: float,
            zones=[],
        ):
        super().__init__("NusseltAndWall", zones)

        self.const = const
        self.coeff = coeff
        self.expRe = expRe
        self.expPr = expPr
        self.addH = addH


    def __repr__(self, depth=0):
        self.__setitem__("const", self.const)
        self.__setitem__("coeff", self.coeff)
        self.__setitem__("expRe", self.expRe)
        self.__setitem__("expPr", self.expPr)
        self.__setitem__("addH", self.addH)

        return super().__repr__(depth)


class NusseltWallAndHfromFMU(NusseltAndWall):
    """
    Nusselt Heat Transfer Model with an additional wall

    Parameters
    ----------
    const : float
    coeff : float
    expRe : float
    expPr : float
    addH : float
    HNameFromFMU : str
        Name of the FMI representing the heat transfer coefficient as a scalar.
    """
    def __init__(
            self,
            const,
            coeff,
            expRe,
            expPr,
            addH,
            HNameFromFMU: str,
            zones=[]
        ):
        super().__init__(const, coeff, expRe, expPr, addH, zones)

        self.type = "NusseltWallAndHfromFMU"
        self.HNameFromFMU = HNameFromFMU


    def __repr__(self, depth=0):
        self.__setitem__("HNameFromFMU", self.HNameFromFMU)

        return super().__repr__(depth)


class Shah(HeatTransferModel):
    """
    Shah Heat Transfer Model

    Implementation Notes
    --------------------

    The original pool boiling heat transfer coefficient by Shah is the following
    form:

    .. math::
        htcPB = C * q^n * pR^m

    with C, n, m being constants, pR being the reduced pressure (i.e. the
    ration of the fluid pressure to the fluid pressure at its critical point)
    and q being the heat flux between wall and fluid. Needless to say, this form
    is not particularly suitable for a numerical implementation in its current
    form as htcPB is required to compute q in the first place. In general:

    .. math::
        q = htc*(Twall-Tf)

    with htc being the total heat transfer coefficient between wall and fluid,
    Twall being the wall temperature and Tf the fluid temperature. By assuming
    that htc ~ htcPB (i.e. that most of the heat flux is due to the pool
    boiling mechanism), one can tha equate:

    .. math::
        htcPB = C * (htcPB*(Twall-Tf))^n * pR^m

    which can be re-arranged to isolate htcPB on the left hand side as:

    .. math::
        htcPB = (C * (Twall-Tf)^n * pR^m)^(1.0/(1.0-n))

    Which is how the correlation is implemented here.
    The assumption of htc ~ htcPB can be rather crude but it is widely employed
    to yiled a feasable numerical implementation in other computer codes too.
    The TRACE code is an example, even though the pool boiling correlation used
    in TRACE is the Gorenflo correlation, not the Shah one. Nonetheless, it is
    still a correlation in which htcPB is a function of the wall heat flux

    A footnote on the parameters C, m, n. Shah provided these values:

    C = 13.7, m = 0.22 if pR < 1e-3
    C = 6.9, m = 0.12 if pR > 1e-3
    n = 0.7 always

    Clearly the discotinuity at pR = 1e-3 should be avoided for numerical
    stabilitiy. For this reason, the values of C and m are linearly interpolated
    in the range (pR0, pR0+deltaPR). Currently, pR0 = 5e-4 and deltaPR = 1e-3.
    The choice of these values is arbitrary and should/could change in the future.


    Parameters
    ----------
    useExplicitHeatFlux : bool
        Default False
    """

    def __init__(
            self,
            useExplicitHeatFlux: bool=False,
            zones=[],
        ):
        super().__init__("Shah", zones)

        self.useExplicitHeatFlux = useExplicitHeatFlux


    def __repr__(self, depth=0):
        self.__setitem__("useExplicitHeatFlux", self.useExplicitHeatFlux)

        return super().__repr__(depth)


class HeatTransferByRegime(HeatTransferModel):
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

        for regimeName, heatTransferModel in self.regimes:
            self.__setitem__(regimeName, heatTransferModel)

        return super().__repr__(depth)

    def add_regime(self, heatTransferModel: HeatTransferModel):
        check_type("heatTransferModel", heatTransferModel, HeatTransferModel)
        regimeName = heatTransferModel.zones[0]
        heatTransferModel.zones = []
        self.regimes.append((regimeName, heatTransferModel))


class NoKazimiHeatTransferModel(HeatTransferModel):
    """
    NoKazimi model for the liquid side of a liquid-vapour heat transfer
    coefficient. Can be only used if a phaseChangeModel is specified as this
    models relies on latent heat.
    """
    def __init__(self, phaseName: str):
        super().__init__("NoKazimi", [phaseName])


class NusseltReynoldsPrandtlPowerFluidFluid(HeatTransferModel):
    """
    Heat transfer coefficient for forced convective flows that is computed
    starting from the Nusselt, which is turn obtained via correlation in the
    form:

        Nu = A_+B_*(Re^C_)*(Pr^D_)

    with Re, Pr being the Reynolds (superficial in two-phase) and Prandlt
    of the fluid, and A, B, C, D selectable consts

    Parameters
    ----------
    const : float
    coeff : float
    expRe : float
    expPr : float
    phaseName : str
    """
    def __init__(
            self,
            const: float,
            coeff: float,
            expRe: float,
            expPr: float,
            phaseName: str,
        ):
        super().__init__("NusseltReynoldsPrandtlPower", [phaseName])

        self.const = const
        self.coeff = coeff
        self.expRe = expRe
        self.expPr = expPr


    def __repr__(self, depth: int=0):
        self.__setitem__("const", self.const)
        self.__setitem__("coeff", self.coeff)
        self.__setitem__("expRe", self.expRe)
        self.__setitem__("expPr", self.expPr)

        return super().__repr__(depth)


class ConstantHeatTransfer(HeatTransferModel):
    def __init__(
            self,
            value: float,
            phaseName: str,
        ):
        super().__init__("constant", [phaseName])

        self.value = value

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value) -> None:
        check_type("value", value, (int, float))
        self._value = value
        self.__setitem__("value", value)


class FlowEnhancementFactor(OpenFOAMDict):
    def __init__(
            self,
            type: str
        ):
        super().__init__()

        self.type = type


    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        # check_value("type", type, _HEAT_TRANSFER_MODEL_TYPES)
        self._type_ = type
        self.__setitem__("type", type)


class CobraTfFlowEnhancementFactor(FlowEnhancementFactor):
    def __init__(self):
        super().__init__("COBRA-TF")


class RezkallahSimsFlowEnhancementFactor(FlowEnhancementFactor):
    def __init__(self, exp: float):
        super().__init__("COBRA-TF")

        self.exp = exp

    @property
    def exp(self):
        return self._exp

    @exp.setter
    def exp(self, exp) -> None:
        check_type("exp", exp, (float, int))
        self._exp = exp
        self.__setitem__("exp", exp)


class SuppressionFactor(OpenFOAMDict):
    def __init__(
            self,
            type: str
        ):
        super().__init__()

        self.type = type

    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        # check_value("type", type, _HEAT_TRANSFER_MODEL_TYPES)
        self._type_ = type
        self.__setitem__("type", type)


class CobraTfSuppressionFactor(SuppressionFactor):
    def __init__(self):
        super().__init__("COBRA-TF")


class ChenSuppressionFactor(SuppressionFactor):
    def __init__(self):
        super().__init__("Chen")


class SuperpositionNucleateBoiling(HeatTransferModel):
    """
    Heat transfer coefficient for both single-phase convective scenarios that
    evolve into nucleate boiling ones, but no sub-cooled boiling region in
    between.
    """
    def __init__(
            self,
            forcedConvectionModel: HeatTransferModel,
            poolBoilingModel: HeatTransferModel,
            flowEnhancementFactorModel: FlowEnhancementFactor,
            suppressionFactorModel: SuppressionFactor,
            zones: list[str]
        ):
        super().__init__("superpositionNucleateBoiling", zones=zones)

        self.forcedConvectionModel = forcedConvectionModel
        self.poolBoilingModel = poolBoilingModel
        self.flowEnhancementFactorModel = flowEnhancementFactorModel
        self.suppressionFactorModel = suppressionFactorModel


    def __repr__(self, depth: int=0):
        self.__setitem__("forcedConvection", self.forcedConvectionModel)
        self.__setitem__("poolBoiling", self.poolBoilingModel)
        self.__setitem__("flowEnhancementFactor", self.flowEnhancementFactorModel)
        self.__setitem__("suppressionFactor", self.suppressionFactorModel)

        return super().__repr__(depth)


    @property
    def forcedConvectionModel(self):
        return self._forcedConvectionModel

    @forcedConvectionModel.setter
    def forcedConvectionModel(self, forcedConvectionModel) -> None:
        check_type("forcedConvectionModel", forcedConvectionModel, HeatTransferModel)
        self._forcedConvectionModel = forcedConvectionModel

    @property
    def poolBoilingModel(self):
        return self._poolBoilingModel

    @poolBoilingModel.setter
    def poolBoilingModel(self, poolBoilingModel) -> None:
        check_type("poolBoilingModel", poolBoilingModel, HeatTransferModel)
        self._poolBoilingModel = poolBoilingModel

    @property
    def flowEnhancementFactorModel(self):
        return self._flowEnhancementFactorModel

    @flowEnhancementFactorModel.setter
    def flowEnhancementFactorModel(self, flowEnhancementFactorModel) -> None:
        check_type("flowEnhancementFactorModel", flowEnhancementFactorModel, FlowEnhancementFactor)
        self._flowEnhancementFactorModel = flowEnhancementFactorModel

    @property
    def suppressionFactorModel(self):
        return self._suppressionFactorModel

    @suppressionFactorModel.setter
    def suppressionFactorModel(self, suppressionFactorModel) -> None:
        check_type("suppressionFactorModel", suppressionFactorModel, SuppressionFactor)
        self._suppressionFactorModel = suppressionFactorModel


class NucleateBoilingOnsetModel(OpenFOAMDict):
    def __init__(
            self,
            type: str
        ):
        super().__init__()

        self.type = type

    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        # check_value("type", type, _HEAT_TRANSFER_MODEL_TYPES)
        self._type_ = type
        self.__setitem__("type", type)


class BasuNucleateBoilingOnsetModel(NucleateBoilingOnsetModel):
    def __init__(
            self,
            surfaceTension: float,
            contactAngle: float
        ):
        super().__init__("Basu")

        self.surfaceTension = surfaceTension
        self.contactAngle = contactAngle

    @property
    def surfaceTension(self):
        return self._surfaceTension

    @surfaceTension.setter
    def surfaceTension(self, surfaceTension) -> None:
        check_type("surfaceTension", surfaceTension, (float, int))
        self._surfaceTension = surfaceTension
        self.__setitem__("surfaceTension", surfaceTension)

    @property
    def contactAngle(self):
        return self._contactAngle

    @contactAngle.setter
    def contactAngle(self, contactAngle) -> None:
        check_type("contactAngle", contactAngle, (float, int))
        self._contactAngle = contactAngle
        self.__setitem__("contactAngle", contactAngle)


class SubCooledBoilingFractionModel(OpenFOAMDict):
    def __init__(
            self,
            type: str
        ):
        super().__init__()

        self.type = type

    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        # check_value("type", type, _HEAT_TRANSFER_MODEL_TYPES)
        self._type_ = type
        self.__setitem__("type", type)


class SahaZuberSubCooledBoilingFractionModel(SubCooledBoilingFractionModel):
    def __init__(self):
        super().__init__("SahaZuber")


class CriticalHeatFluxModel(OpenFOAMDict):
    def __init__(
            self,
            type: str
        ):
        super().__init__()

        self.type = type

    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        # check_value("type", type, _HEAT_TRANSFER_MODEL_TYPES)
        self._type_ = type
        self.__setitem__("type", type)


class ConstantCHFCriticalHeatFluxModel(CriticalHeatFluxModel):
    def __init__(
            self,
            value: float
        ):
        super().__init__("constantCHF")

        self.value = value

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value) -> None:
        check_type("value", value, (float, int))
        self._value = value
        self.__setitem__("value", value)


class LeidenfrostModel(OpenFOAMDict):
    def __init__(
            self,
            type: str
        ):
        super().__init__()

        self.type = type

    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        # check_value("type", type, _HEAT_TRANSFER_MODEL_TYPES)
        self._type_ = type
        self.__setitem__("type", type)


class GroeneveldStewartLeidenfrostModel(LeidenfrostModel):
    def __init__(
            self,
            criticalPressure: float
        ):
        super().__init__("GroeneveldStewart")

        self.criticalPressure = criticalPressure

    @property
    def criticalPressure(self):
        return self._criticalPressure

    @criticalPressure .setter
    def criticalPressure(self, criticalPressure) -> None:
        check_type("criticalPressure", criticalPressure, (float, int))
        self._criticalPressure = criticalPressure
        self.__setitem__("criticalPressure", criticalPressure)


class AnnularFlowModel(OpenFOAMDict):
    def __init__(
            self,
            type: str
        ):
        super().__init__()

        self.type = type

    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        # check_value("type", type, _HEAT_TRANSFER_MODEL_TYPES)
        self._type_ = type
        self.__setitem__("type", type)


class CachardLiquidAnnularFlowModel(AnnularFlowModel):
    def __init__(
            self,
            wallEmissivity: float,
            liquidEmissivity: float
        ):
        super().__init__("CachardLiquid")

        self.wallEmissivity = wallEmissivity
        self.liquidEmissivity = liquidEmissivity

    @property
    def wallEmissivity(self):
        return self._wallEmissivity

    @wallEmissivity .setter
    def wallEmissivity(self, wallEmissivity) -> None:
        check_type("wallEmissivity", wallEmissivity, (float, int))
        self._wallEmissivity = wallEmissivity
        self.__setitem__("wallEmissivity", wallEmissivity)

    @property
    def liquidEmissivity(self):
        return self._liquidEmissivity

    @liquidEmissivity .setter
    def liquidEmissivity(self, liquidEmissivity) -> None:
        check_type("liquidEmissivity", liquidEmissivity, (float, int))
        self._liquidEmissivity = liquidEmissivity
        self.__setitem__("liquidEmissivity", liquidEmissivity)


class CachardVapourAnnularFlowModel(AnnularFlowModel):
    def __init__(self):
        super().__init__("CachardVapour")


class MultiRegimeBoilingTRACECHF(HeatTransferModel):
    """
    Heat transfer coefficient for multi-regime boiling from TRACE CHF.
    """
    def __init__(
            self,
            forcedConvectionModel: HeatTransferModel,
            poolBoilingModel: HeatTransferModel,
            flowEnhancementFactorModel: FlowEnhancementFactor,
            suppressionFactorModel: SuppressionFactor,
            nucleateBoilingOnsetModel: NucleateBoilingOnsetModel,
            subCooledBoilingFractionModel: SubCooledBoilingFractionModel,
            criticalHeatFluxModel: CriticalHeatFluxModel,
            leidenfrostModel: LeidenfrostModel,
            annularFlowModel: AnnularFlowModel,
            zones: list[str]
        ):
        super().__init__("multiRegimeBoilingTRACECHF", zones=zones)

        self.forcedConvectionModel = forcedConvectionModel
        self.poolBoilingModel = poolBoilingModel
        self.flowEnhancementFactorModel = flowEnhancementFactorModel
        self.suppressionFactorModel = suppressionFactorModel
        self.nucleateBoilingOnsetModel = nucleateBoilingOnsetModel
        self.subCooledBoilingFractionModel = subCooledBoilingFractionModel
        self.criticalHeatFluxModel = criticalHeatFluxModel
        self.leidenfrostModel = leidenfrostModel
        self.annularFlowModel = annularFlowModel


    def __repr__(self, depth: int=0):
        self.__setitem__("forcedConvectionModel", self.forcedConvectionModel)
        self.__setitem__("poolBoilingModel", self.poolBoilingModel)
        self.__setitem__("flowEnhancementFactorModel", self.flowEnhancementFactorModel)
        self.__setitem__("suppressionFactorModel", self.suppressionFactorModel)
        self.__setitem__("nucleateBoilingOnsetModel", self.nucleateBoilingOnsetModel)
        self.__setitem__("subCooledBoilingFractionModel", self.subCooledBoilingFractionModel)
        self.__setitem__("criticalHeatFluxModel", self.criticalHeatFluxModel)
        self.__setitem__("leidenfrostModel", self.leidenfrostModel)
        self.__setitem__("annularFlowModel", self.annularFlowModel)

        return super().__repr__(depth)


    @property
    def forcedConvectionnModel(self):
        return self._forcedConvectionnModel

    @forcedConvectionnModel.setter
    def forcedConvectionnModel(self, forcedConvectionnModel) -> None:
        check_type("forcedConvectionnModel", forcedConvectionnModel, HeatTransferModel)
        self._forcedConvectionnModel = forcedConvectionnModel

    @property
    def poolBoilingnModel(self):
        return self._poolBoilingnModel

    @poolBoilingnModel.setter
    def poolBoilingnModel(self, poolBoilingnModel) -> None:
        check_type("poolBoilingnModel", poolBoilingnModel, HeatTransferModel)
        self._poolBoilingnModel = poolBoilingnModel

    @property
    def flowEnhancementFactorModel(self):
        return self._flowEnhancementFactorModel

    @flowEnhancementFactorModel.setter
    def flowEnhancementFactorModel(self, flowEnhancementFactorModel) -> None:
        check_type("flowEnhancementFactorModel", flowEnhancementFactorModel, FlowEnhancementFactor)
        self._flowEnhancementFactorModel = flowEnhancementFactorModel

    @property
    def suppressionFactorModel(self):
        return self._suppressionFactorModel

    @suppressionFactorModel.setter
    def suppressionFactorModel(self, suppressionFactorModel) -> None:
        check_type("suppressionFactorModel", suppressionFactorModel, SuppressionFactor)
        self._suppressionFactorModel = suppressionFactorModel

    @property
    def nucleateBoilingOnsetModel(self):
        return self._nucleateBoilingOnsetModel

    @nucleateBoilingOnsetModel.setter
    def nucleateBoilingOnsetModel(self, nucleateBoilingOnsetModel) -> None:
        check_type("nucleateBoilingOnsetModel", nucleateBoilingOnsetModel, NucleateBoilingOnsetModel)
        self._nucleateBoilingOnsetModel = nucleateBoilingOnsetModel

    @property
    def subCooledBoilingFractionModel(self):
        return self._subCooledBoilingFractionModel

    @subCooledBoilingFractionModel.setter
    def subCooledBoilingFractionModel(self, subCooledBoilingFractionModel) -> None:
        check_type("subCooledBoilingFractionModel", subCooledBoilingFractionModel, SubCooledBoilingFractionModel)
        self._subCooledBoilingFractionModel = subCooledBoilingFractionModel

    @property
    def criticalHeatFluxModel(self):
        return self._criticalHeatFluxModel

    @criticalHeatFluxModel.setter
    def criticalHeatFluxModel(self, criticalHeatFluxModel) -> None:
        check_type("criticalHeatFluxModel", criticalHeatFluxModel, CriticalHeatFluxModel)
        self._criticalHeatFluxModel = criticalHeatFluxModel

    @property
    def leidenfrostModel(self):
        return self._leidenfrostModel

    @leidenfrostModel.setter
    def leidenfrostModel(self, leidenfrostModel) -> None:
        check_type("leidenfrostModel", leidenfrostModel, LeidenfrostModel)
        self._leidenfrostModel = leidenfrostModel

    @property
    def annularFlowModel(self):
        return self._annularFlowModel

    @annularFlowModel.setter
    def annularFlowModel(self, annularFlowModel) -> None:
        check_type("annularFlowModel", annularFlowModel, AnnularFlowModel)
        self._annularFlowModel = annularFlowModel


class MultiRegimeBoilingVapourTRACE(HeatTransferModel):
    """
    Heat transfer coefficient for multi-regime boiling for vapour from TRACE.
    """
    def __init__(
            self,
            annularFlowModel: AnnularFlowModel,
            zones: list[str]
        ):
        super().__init__("multiRegimeBoilingVapourTRACE", zones=zones)

        self.annularFlowModel = annularFlowModel


    def __repr__(self, depth: int=0):
        self.__setitem__("annularFlowModel", self.annularFlowModel)

        return super().__repr__(depth)

    @property
    def annularFlowModel(self):
        return self._annularFlowModel

    @annularFlowModel.setter
    def annularFlowModel(self, annularFlowModel) -> None:
        check_type("annularFlowModel", annularFlowModel, AnnularFlowModel)
        self._annularFlowModel = annularFlowModel
