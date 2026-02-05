# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import offbeat_define
from foamForNuclear.common import OffbeatDict

@offbeat_define
class NeutronicsSubSolver(OffbeatDict):
    """
    Base neutronics solver class. It neglects entirely the evolution of the neutron
    flux over time.

    No additional neutron-related field is created and no neutronics equation
    is solved.
    """
    TYPE: ClassVar[str] = 'none'

@offbeat_define
class Diffusion(NeutronicsSubSolver):
    """
    Neutronics solver class that obtains the neutron flux distribution from the
    solution of a one-group neutron diffusion equation using fitted parameters.
    """
    TYPE: ClassVar[str] = 'diffusion'
