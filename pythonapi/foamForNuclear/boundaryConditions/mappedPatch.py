from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import Vector


_SAMPLE_MODE_TYPES = {
    "nearestCell", "nearestOnlyCell", "nearestPatchFace", "nearestPatchFaceAMI",
    "nearestFace", "nearestPatchPoint"
}
_OFFSET_MODE_TYPES = {"uniform", "nonuniform", "normal"}

class MappedPatch(Patch):
    """
    The mappedField is a generic boundary condition that provides a
    self-contained version of the mapped condition. It does not use information
    on the patch; instead it holds the data locally.

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
        - `nearestPatchFaceAMI`: nearest face on selected patch (patches need
            not conform; uses AMI interpolation)
        - `nearestFace`: nearest boundary face on any patch
        - `nearestPatchPoint`: nearest patch point (for coupled points this
            might be any of the points so you have to guarantee the point data
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
            type: str,
            value: float,
            sampleRegion: str,
            sampleMode: str,
            samplePatch: str,
            offsetMode: str=None,
            offset: Vector=None,
            coupleGroup: str=None,
            field: str=None,
            setAverage: bool=None,
            average: float=None
        ):
        super().__init__(type=type, value=value)

        self.sampleRegion = sampleRegion
        self.sampleMode = sampleMode
        self.samplePatch = samplePatch
        self.coupleGroup = coupleGroup
        self.offsetMode = offsetMode
        self.offset = offset
        self.field = field
        self.setAverage = setAverage
        self.average = average


    @property
    def sampleRegion(self):
        return self._sampleRegion

    @sampleRegion.setter
    def sampleRegion(self, sampleRegion) -> None:
        check_type("sampleRegion", sampleRegion, str)
        self._sampleRegion = sampleRegion
        self.__setitem__('sampleRegion', sampleRegion)

    @property
    def sampleMode(self):
        return self._sampleMode

    @sampleMode.setter
    def sampleMode(self, sampleMode) -> None:
        check_type("sampleMode", sampleMode, str)
        check_value("sampleMode", sampleMode, _SAMPLE_MODE_TYPES)
        self._sampleMode = sampleMode
        self.__setitem__('sampleMode', sampleMode)

    @property
    def samplePatch(self):
        return self._samplePatch

    @samplePatch.setter
    def samplePatch(self, samplePatch) -> None:
        check_type("samplePatch", samplePatch, str)
        self._samplePatch = samplePatch
        self.__setitem__('samplePatch', samplePatch)

    @property
    def coupleGroup(self):
        return self._coupleGroup

    @coupleGroup.setter
    def coupleGroup(self, coupleGroup) -> None:
        check_type("coupleGroup", coupleGroup, str, none_ok=True)
        self._coupleGroup = coupleGroup
        if (coupleGroup is not None):
            self.__setitem__('coupleGroup', coupleGroup)

    @property
    def offsetMode(self):
        return self._offsetMode

    @offsetMode.setter
    def offsetMode(self, offsetMode) -> None:
        check_type("offsetMode", offsetMode, str, none_ok=True)
        if (offsetMode is not None):
            check_value("offsetMode", offsetMode.split()[0], _OFFSET_MODE_TYPES)
            self.__setitem__('offsetMode', offsetMode)
        self._offsetMode = offsetMode

    @property
    def offset(self):
        return self._offset

    @offset.setter
    def offset(self, offset) -> None:
        check_type("offset", offset, Vector, none_ok=True)
        self._offset = offset
        if (offset is not None):
            self.__setitem__('offset', offset)

    @property
    def field(self):
        return self._field

    @field.setter
    def field(self, field) -> None:
        check_type("field", field, str, none_ok=True)
        self._field = field
        if (field is not None):
            self.__setitem__('field', field)

    @property
    def setAverage(self):
        return self._setAverage

    @setAverage.setter
    def setAverage(self, setAverage) -> None:
        check_type("setAverage", setAverage, bool, none_ok=True)
        self._setAverage = setAverage
        if (setAverage is not None):
            self.__setitem__('setAverage', setAverage)

    @property
    def average(self):
        return self._average

    @average.setter
    def average(self, average) -> None:
        check_type("average", average, (float, int), none_ok=True)
        self._average = average
        if (average is not None):
            self.__setitem__('average', average)
