from collections.abc import Mapping

from foamForNuclear.common import FoamForNuclearDict, List, OpenFOAMDict
from foamForNuclear.checkvalue import check_type
from foamForNuclear._attrs_tools import ffn_define, _to_List_str
from typing import ClassVar
from attrs import field, validators as v


from .annular_flow import AnnularFlowModel
from .chf import CriticalHeatFluxModel
from .enhancement import FlowEnhancementFactor
from .leidenfrost import LeidenfrostModel
from .onb import NucleateBoilingOnsetModel
from .subcooled_fraction import SubCooledBoilingFractionModel
from .suppression import SuppressionFactorModel


@ffn_define
class HeatTransferModel(FoamForNuclearDict):
    """
    Base class for heat-transfer.
    """
    TYPE: ClassVar[str] = "none"
    zones: list | List = field(factory=list, metadata={"ffn_internal": True}, converter=_to_List_str)

    def __repr__(self, depth = 0):
        textZones = ""
        if (len(self.zones) > 0):
            textZones = "\"" + ':'.join(self.zones) + "\""

        return textZones + super().__repr__(depth)


@ffn_define
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
    TYPE: ClassVar[str] = "Gorenflo"
    absoluteSurfaceRoughness: int | float = 4e-7
    useExplicitHeatFlux: bool = False


@ffn_define
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
    TYPE: ClassVar[str] = "NusseltAndWall"
    const: int | float
    coeff: int | float
    expRe: int | float
    expPr: int | float
    addH:  int | float


@ffn_define
class NusseltWallAndHfromFMU(HeatTransferModel):
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
    TYPE: ClassVar[str] = "NusseltWallAndHfromFMU"
    const: int | float
    coeff: int | float
    expRe: int | float
    expPr: int | float
    addH: int | float
    HNameFromFMU: str


@ffn_define
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
    TYPE: ClassVar[str] = "Shah"
    useExplicitHeatFlux: bool = False


@ffn_define
class ByRegime(HeatTransferModel):
    TYPE: ClassVar[str] = "byRegime"
    regimeMap: str
    regimes: list[HeatTransferModel] = field(
        factory=list,
        metadata={"ffn_internal": True},
    )

    def _update_fields(self) -> None:
        super()._update_fields()

        for reg in self.regimes:
            if reg is None:
                continue

            # --- Case 1: flat dict API ---
            if isinstance(reg, Mapping):
                name = reg.get("name")
                if not name:
                    raise ValueError("Each regime dict must define a 'name'")

                # everything except 'name' goes into the OpenFOAM sub-dict
                items = {k: v for k, v in reg.items() if k != "name"}
                if not items:
                    raise ValueError(f"Regime '{name}' has no entries")

                reg_of = OpenFOAMDict(name=name, items=items)

            # --- Case 2: OpenFOAMDict-like object ---
            else:
                name = getattr(reg, "name", None)
                if not name:
                    raise ValueError("Each regime must have a 'name'")
                reg_of = reg

            self._of[name] = reg_of

    # def __repr__(self, depth=0):
    #     self.__setitem__("regimeMap", self.regimeMap)

    #     for regimeName, heatTransferModel in self.regimes:
    #         self.__setitem__(regimeName, heatTransferModel)

    #     return super().__repr__(depth)

    def add_regime(self, heatTransferModel: HeatTransferModel):
        check_type("heatTransferModel", heatTransferModel, HeatTransferModel)
        regimeName = heatTransferModel.zones[0]
        heatTransferModel.zones = []
        self.regimes.append((regimeName, heatTransferModel))


@ffn_define
class NoKazimiHeatTransferModel(HeatTransferModel):
    """
    NoKazimi model for the liquid side of a liquid-vapour heat transfer
    coefficient. Can be only used if a phaseChangeModel is specified as this
    models relies on latent heat.
    """
    TYPE: ClassVar[str] = "NoKazimi"


@ffn_define
class NusseltReynoldsPrandtlPower(HeatTransferModel):
    """
    Heat transfer coefficient for forced convective flows that is computed
    starting from the Nusselt, which is turn obtained via correlation in the
    form:

        Nu = A_+B_*(Re^C_)*(Pr^D_)

    with Re, Pr being the Reynolds (superficial in two-phase) and Prandlt
    of the fluid, and A, B, C, D selectable consts

    Parameters
    ----------
    const : int | float
    coeff : int | float
    expRe : int | float
    expPr : int | float
    expTc : Optional[int | float]
    phaseName : Optional[str] # Only for fluid-fluid interaction
    """
    TYPE: ClassVar[str] = "NusseltReynoldsPrandtlPower"
    const: int | float
    coeff: int | float
    expRe: int | float
    expPr: int | float
    expTc: int | float | None = None
    phaseName: str | None = field(default=None, metadata={"ffn_internal": True})


@ffn_define
class Constant(HeatTransferModel):
    TYPE: ClassVar[str] = "constant"
    value: int | float


@ffn_define
class SuperpositionNucleateBoiling(HeatTransferModel):
    """
    Heat transfer coefficient for both single-phase convective scenarios that
    evolve into nucleate boiling ones, but no sub-cooled boiling region in
    between.
    """
    TYPE: ClassVar[str] = "superpositionNucleateBoiling"
    forcedConvection: HeatTransferModel = field(factory=HeatTransferModel)
    poolBoiling: HeatTransferModel = field(factory=HeatTransferModel)
    flowEnhancementFactor: FlowEnhancementFactor = field(factory=FlowEnhancementFactor)
    suppressionFactor: SuppressionFactorModel = field(factory=SuppressionFactorModel)


@ffn_define
class MultiRegimeBoilingTRACECHF(HeatTransferModel):
    """
    Heat transfer coefficient for multi-regime boiling from TRACE CHF.
    """
    TYPE: ClassVar[str] = "multiRegimeBoilingTRACECHF"
    forcedConvectionModel: HeatTransferModel = field(factory=HeatTransferModel)
    poolBoilingModel: HeatTransferModel = field(factory=HeatTransferModel)
    flowEnhancementFactorModel: FlowEnhancementFactor = field(factory=FlowEnhancementFactor)
    suppressionFactorModel: SuppressionFactorModel = field(factory=SuppressionFactorModel)
    nucleateBoilingOnsetModel: NucleateBoilingOnsetModel = field(factory=NucleateBoilingOnsetModel)
    subCooledBoilingFractionModel: SubCooledBoilingFractionModel = field(factory=SubCooledBoilingFractionModel)
    criticalHeatFluxModel: CriticalHeatFluxModel = field(factory=CriticalHeatFluxModel)
    leidenfrostModel: LeidenfrostModel = field(factory=LeidenfrostModel)
    annularFlowModel: AnnularFlowModel = field(factory=AnnularFlowModel)


@ffn_define
class MultiRegimeBoilingVapourTRACE(HeatTransferModel):
    """
    Heat transfer coefficient for multi-regime boiling for vapour from TRACE.
    """
    TYPE: ClassVar[str] = "multiRegimeBoilingVapourTRACE"
    annularFlowModel: AnnularFlowModel = field(factory=AnnularFlowModel)