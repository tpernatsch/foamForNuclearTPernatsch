from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type


class AlbedoSP3(Patch):
    """
    Albedo boundary condition for SP3 or diffusion calculations. Please note
    that the boundary condition needs to be set both for first and second
    moments in SP3.

    Parameters
    ----------
    value : float | int
        Value at the patch set uniformly.
    alpha : float
        Ratio between outcoming and incoming flux.
    gamma : float
        Defined as `(1-alpha)/(1+alpha)/2`.
    diffCoeffName : str
        Name of the diffusion coefficient field (default `Dalbedo`).
    fluxStarAlbedo : str
        Name of the other moment (default `fluxStarAlbedo`).
    forSecondMoment : bool
        True in case it is a condition for a second moment flux (for SP3
        calculations) (default `False`).
    """
    def __init__(
            self,
            value: float,
            alpha: float=None,
            gamma: float=None,
            forSecondMoment: bool=False,
            diffCoeffName: str="Dalbedo",
            fluxStarAlbedo: str="fluxStarAlbedo",
        ):
        super().__init__(type="albedoSP3", value=value)

        self.gamma = gamma
        self.alpha = alpha
        if (self.alpha is not None):
            self.gamma = self.getGamma(alpha)
        elif (self.gamma is not None):
            self.alpha = self.getAlpha(gamma)

        self.diffCoeffName = diffCoeffName
        self.fluxStarAlbedo = fluxStarAlbedo
        self.forSecondMoment = forSecondMoment


    def __repr__(self, depth = 0):
        self.__setitem__('gamma', self.gamma)
        self.__setitem__('diffCoeffName', self.diffCoeffName)
        self.__setitem__('fluxStarAlbedo', self.fluxStarAlbedo)
        self.__setitem__('forSecondMoment', self.forSecondMoment)

        return super().__repr__(depth)


    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, alpha) -> None:
        check_type("alpha", alpha, (float, int), none_ok=True)
        self._alpha = alpha
        if (alpha is not None):
            self._gamma = self.getGamma(alpha)


    @property
    def gamma(self):
        return self._gamma

    @gamma.setter
    def gamma(self, gamma) -> None:
        check_type("gamma", gamma, (float, int), none_ok=True)
        self._gamma = gamma
        if (gamma is not None):
            self._alpha = self.getAlpha(gamma)


    @property
    def diffCoeffName(self):
        return self._diffCoeffName

    @diffCoeffName.setter
    def diffCoeffName(self, diffCoeffName) -> None:
        check_type("diffCoeffName", diffCoeffName, str)
        self._diffCoeffName = diffCoeffName


    @property
    def fluxStarAlbedo(self):
        return self._fluxStarAlbedo

    @fluxStarAlbedo.setter
    def fluxStarAlbedo(self, fluxStarAlbedo) -> None:
        check_type("fluxStarAlbedo", fluxStarAlbedo, str)
        self._fluxStarAlbedo = fluxStarAlbedo


    @property
    def forSecondMoment(self):
        return self._forSecondMoment

    @forSecondMoment.setter
    def forSecondMoment(self, forSecondMoment) -> None:
        check_type("forSecondMoment", forSecondMoment, bool)
        self._forSecondMoment = forSecondMoment


    def getAlpha(self, gamma):
        return((1 - 2*gamma) / (1 + 2*gamma))

    def getGamma(self, alpha):
        return(0.5 * (1 - alpha) / (1 + alpha))
