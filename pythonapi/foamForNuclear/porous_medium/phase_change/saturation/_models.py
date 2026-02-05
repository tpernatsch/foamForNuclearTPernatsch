import numpy as np
from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import OffbeatDict, List, OpenFOAMDict
from foamForNuclear._attrs_tools import offbeat_define, _to_List_str
from typing import ClassVar
from attrs import field, validators as v



@offbeat_define
class SaturationModel(OffbeatDict):
    TYPE: ClassVar[str] = "none"


@offbeat_define
class BrowningPotter(SaturationModel):
    """
    Saturation model based on results by Browning and Potter, refer to
    https://www.ne.anl.gov/eda/ANL-RE-95-2.pdf. The saturation pressure form is
    lnPSat = A + B/TSat + C*ln(TSat), which is not invertible in TSat. TSat is
    inverted by approximateding ln(TSat) by a second order polyinomial.
    """
    TYPE: ClassVar[str] = "BrowningPotter"


@offbeat_define
class ConstantTemperature(SaturationModel):
    """
    Constant saturation pressure and temperature.

    Parameters
    ----------
    value : float
        Saturation temperature
    """
    TYPE: ClassVar[str] = "constantTemperature"


@offbeat_define
class Water(SaturationModel):
    """
    Saturation model for water in the 0.01-350 C range based on an
    interpolation of NIST data in the form

    .. math::
        Psat = A (T-B)^C

    The NIST data can be found at this `link
    <https://www.nist.gov/system/files/documents/srd/NISTIR5078-Tab1.pdf>`_.
    """
    TYPE: ClassVar[str] = "water"


@offbeat_define
class WaterTRACE(SaturationModel):
    """
    Saturation model for water in the 0.01-350 C range based on `TRACE
    <https://www.nrc.gov/docs/ML1200/ML120060218.pdf>`_.
    """
    TYPE: ClassVar[str] = "waterTRACE"

