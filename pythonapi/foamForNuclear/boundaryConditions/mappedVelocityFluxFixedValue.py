from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.boundaryConditions.mappedPatch import MappedPatch
from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import Vector



class MappedVelocityFluxFixedValue(MappedPatch):
    """
    The mappedVelocityFluxFixedValue is a velocity boundary condition that maps
    the velocity and flux from a neighbour patch to this patch.

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
    nbrPhi : str
        Name of the neighbour flux setting the inlet mass flux (default `phi`).
    phi : str
        Name of the local mass flux (default `phi`).
    rho : str
        Name of the density field used to normalize the mass flux
        (default `rho`).
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
            nbrPhi: str='phi',
            phi: str='phi',
            rho: str='rho',
            offsetMode: str='uniform',
            offset: Vector=Vector(0, 0, 0),
            coupleGroup: str=None,
            field: str=None,
            setAverage: bool=False,
            average: float=0
        ):
        super().__init__(
            type="mappedVelocityFlux",
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

        self.nbrPhi = nbrPhi
        self.phi = phi
        self.rho = rho

    @property
    def nbrPhi(self):
        return self._nbrPhi

    @nbrPhi.setter
    def nbrPhi(self, nbrPhi) -> None:
        check_type("nbrPhi", nbrPhi, str)
        self._nbrPhi = nbrPhi
        self.__setitem__('nbrPhi', nbrPhi)

    @property
    def phi(self):
        return self._phi

    @phi.setter
    def phi(self, phi) -> None:
        check_type("phi", phi, str)
        self._phi = phi
        self.__setitem__('phi', phi)

    @property
    def rho(self):
        return self._rho

    @rho.setter
    def rho(self, rho) -> None:
        check_type("rho", rho, str)
        self._rho = rho
        self.__setitem__('rho', rho)
