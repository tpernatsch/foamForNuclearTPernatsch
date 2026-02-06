import numpy as np
from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import FoamForNuclearDict, List, OpenFOAMDict
from foamForNuclear._attrs_tools import ffn_define, _to_List_str
from typing import ClassVar
from attrs import field, validators as v


@ffn_define
class LatentHeatModel(FoamForNuclearDict):
    TYPE: ClassVar[str] = "none"
    adjust: bool = False


@ffn_define
class FinkLeibowitz(LatentHeatModel):
    """
    This class describes the latent heat of vaporization of liquid sodium
    according to the data presented in the report by Fink & Leibowitz
    (ANL RE 95/2), https://www.ne.anl.gov/eda/ANL-RE-95-2.pdf
    """
    TYPE: ClassVar[str] = "FinkLeibowitz"


@ffn_define
class FromThermophysicalProperties(LatentHeatModel):
    """
    Class that computes latent heat based on the specified enthalpies of
    formation of the fluids. These are specified under the Hf keyword
    in the thermophysical properties of the two fluids. The latent heat is thus
    computed as Hf.vapour - Hf.liquid.
    """
    TYPE: ClassVar[str] = "fromThermophysicalProperties"


@ffn_define
class Water(LatentHeatModel):
    """
    This class describes the latent heat of vaporization of water in the
    0.01-350 C range. It was obtained from a fit (in the form of
    L = A+B*log(C-T) with T being the liquid temperature (K)) of NIST water
    data at saturation found at:
    https://www.nist.gov/system/files/documents/srd/NISTIR5078-Tab1.pdf

    It performs reasonably well up to ~ 365 C as well, considering that the
    critical water temperature is 375.16 C (at which point the latent heat
    becomes 0, but is limited for obvious numerical reasons as this code is not
    meant for sub-supercritical fluid transition simulations)
    """
    TYPE: ClassVar[str] = "water"

