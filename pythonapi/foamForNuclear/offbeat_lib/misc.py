from __future__ import annotations

import attrs as attr
from attrs import field, validators as v
from foamForNuclear._attrs_tools import offbeat_define, _to_List_float, _to_List_str, _is_number_when_flag
from foamForNuclear._attrs_tools import _to_List_tuple_float, _to_List_tuple2_float, _resolve_type, _to_Vector
from typing import ClassVar

import foamForNuclear
from foamForNuclear.common import *
import foamForNuclear.common

#==============================================================================*

@offbeat_define
class ThermoMechanicsCouplingOptions(OffbeatDict):
    correctTFromTH: bool = False
    correctDispForNeutro: bool = False


@offbeat_define
class GlobalOptions(OffbeatDict):
    pinDirection: list[float | int] | None = None
    reactorType: str | None = None
    angularFraction: float | None = None
    linkedFuel: bool = False


#==============================================================================*
# Element transport solver

@offbeat_define
class ElementTransportSolver(OffbeatDict):
    TYPE: ClassVar[str] = "fromLatestTime"


@offbeat_define
class ByListElementTransportSolver(ElementTransportSolver):
    TYPE: ClassVar[str] = "byList"
    solvers: list | List = field(factory=list, converter=_to_List_str)
  

#==============================================================================*
# Stress Analysis

@offbeat_define
class StressAnalysis(OffbeatDict):
    """
    StressAnalysis object

    Parameters
    ----------
    nCorrectors : int | dict[str, int]
    maxOuterIter : int
    referencePairs : List | list | np.ndarray
    D : Vector
    T : float
    relD : float
    relT : float
    """
    
    nCorrectors: dict[str, int] | int = 1
    maxOuterIter: int = 100
    referencePairs: list | List = field(factory=list, converter=_to_List_tuple2_float)
    D: list | List = field(default=attr.Factory(lambda self: List([1e-5, 1e-5, 1e-5]), takes_self=True), converter=_to_List_float)
    T: float = 1e-5
    neutronFlux0: float = 1e-5
    useRelResD: bool = False
    useRelResT: bool = False
    useRelResNeutronFlux0: bool = False
    relD: float | None = field(default=None, validator=_is_number_when_flag("useRelResD"))
    relT: float | None = field(default=None, validator=_is_number_when_flag("useRelResT"))
    relneutronFlux0: float | None = field(default=None, validator=_is_number_when_flag("useRelResNeutronFlux0"))
    absErrD: float | int | None = None
    absErrT: float | int | None = None