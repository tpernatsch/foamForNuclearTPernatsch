# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import offbeat_define
from foamForNuclear.common import OffbeatDict

@offbeat_define
class SliceMapper(OffbeatDict):
    """
    Base sliceMapper class. It defines the strategy used to map mesh cells or regions
    into axial or logical slices.

    Slice mappers are used to aggregate field quantities (e.g. temperature,
    power, burnup) over geometrically or logically defined virtual slices.

    The base class provides the common interface and infrastructure, but if selected
    it deactivate the creation of virtual axial slices, i.e. no
    virtual slices are defined.
    """
    TYPE: ClassVar[str] = 'none'

@offbeat_define
class AutoAxialSlices(SliceMapper):
    """
    Slice mapper class where slices are automaticaly created for each material.


    Options
    -------
    precision : scalar
        Precision used to group cell-center axial coordinates when building slices.
        The axial coordinate is rounded to multiples of `precision`.
        (default: 1e-06; required: False)
    """
    TYPE: ClassVar[str] = 'autoAxialSlices'
    precision: float | int = 1e-06

@offbeat_define
class ByMaterial(SliceMapper):
    """
    Slice mapper class where slices are created material by material.
    The number and height of the slices are provided by the user in the materials
    subdictionary (by default the slice number in a material is 1).
    """
    TYPE: ClassVar[str] = 'byMaterial'

@offbeat_define
class ByPellets(SliceMapper):
    """
    Slice mapper class where slices are created with a similar algorithm as in the
    `byMaterial` mapper class. For fuel materials (e.g. UO2) the slice height
    coincide with the pellet size, introduced by the user in the material subdictionary.
    For non-fuel materials, all cells are agglomerated into 1 slice.

    The algorithm tries to create one slice per pellet, which for 2D or 3D discrete
    pellet simulations is probably the most logical approximation for phenomena such
    as relocation or cracking (e.g. no variation of number of cracks along a pellet)

    For coarser simulations (e.g. in 1D or 2D smeared), the cells will be taller
    than the pellet size. The slice division algorithm will adapt following the
    mesh axial discretization. For extruded meshes, this will create one slice
    per axial layer of cells.
    """
    TYPE: ClassVar[str] = 'byPellets'
