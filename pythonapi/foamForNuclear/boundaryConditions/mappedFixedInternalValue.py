from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.boundaryConditions.mappedPatch import MappedPatch
from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import Vector



class MappedFixedInternalValue(MappedPatch):
    """
    This boundary condition maps the boundary and internal values of a
    neighbour patch field to the boundary and internal values of this.

    Parameters
    ----------
    value : float | int | Vector
        Value at the patch set uniformly
    sampleRegion : str
        Name of the region to sample
    sampleMode : str
        - `nearestCell`: sample cell containing point
        - `nearestOnlyCell`: nearest sample cell (even if not containing point)
        - `nearestPatchFace`: nearest face on selected patch
        - `nearestPatchFaceAMI`: nearest face on selected patch (patches need \
            not conform; uses AMI interpolation)
        - `nearestFace`: nearest boundary face on any patch
        - `nearestPatchPoint`: nearest patch point (for coupled points this \
            might be any of the points so you have to guarantee the point data \
            is synchronised beforehand)

    samplePatch  : str
        If sampleMode is `nearestPatchFace`: patch to find faces of
    offsetMode : str
        Default `uniform`.
        How to supply offset (w.r.t. my patch face centres):

        - `uniform`: single offset vector
        - `nonuniform`: per-face offset vector
        - `normal`: using supplied distance and face normal

    offset : Vector
        Default `Vector(0, 0, 0)`.
        According to `offsetMode` (see above) supply one of offset, offsets or
        distance
    coupleGroup : str
        If sampleMode is `nearestPatchFace`: specify patchgroup to find
        `samplePatch` and `sampleRegion` (if not provided)
    field : str
        Default `None`
    setAverage : bool
        Default `False`
    average : float
        Default `0`
    """
    def __init__(
            self,
            value: float,
            sampleRegion: str,
            sampleMode: str,
            samplePatch: str,
            offsetMode: str='uniform',
            offset: Vector=Vector(0, 0, 0),
            coupleGroup: str=None,
            field: str=None,
            setAverage: bool=False,
            average: float=0
        ):
        super().__init__(
            type="mappedFixedInternalValue",
            value=value,
            sampleRegion=sampleRegion,
            sampleMode=sampleMode,
            samplePatch=samplePatch,
            offsetMode=offsetMode,
            offset=offset,
            coupleGroup=coupleGroup,
            field=field,
            setAverage=setAverage,
            average=average
        )
