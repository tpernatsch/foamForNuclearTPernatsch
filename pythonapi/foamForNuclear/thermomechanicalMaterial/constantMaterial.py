from foamForNuclear.checkvalue import check_type
from foamForNuclear.field import ReducedDimension
from .thermomechanicalMaterial import BaseThermomechanicalMaterial, ConstantMechanicalLaw


class ConstantMaterial(BaseThermomechanicalMaterial):
    def __init__(
            self,
            name: str,
            rho: float,
            Cp: float,
            k: float,
            emissivity: float,
            E: float,
            nu: float,
            G: float,
            alpha: float,
            Tref: float,
            alphaFuel: float=None,
            TFuelRef: float=None,
            alphaCR: float=None,
            TCRRef: float=None,
            rheologyModel: str="elasticity"
        ):
        super().__init__(name=name, material="constant")

        self.rho = rho
        self.Cp = Cp
        self.k = k
        self.emissivity = emissivity
        self.E = E
        self.nu = nu
        self.G = G
        self.alpha = alpha
        self.Tref = Tref
        self.alphaFuel = alphaFuel
        self.TFuelRef = TFuelRef
        self.alphaCR = alphaCR
        self.TCRRef = TCRRef
        self.rheologyModel = rheologyModel


    def __repr__(self, depth = 0):
        self.__setitem__('rho', ConstantMechanicalLaw("rho", ReducedDimension('rho'), self.rho))
        self.__setitem__('Cp', ConstantMechanicalLaw("Cp", ReducedDimension('Cp'), self.Cp))
        self.__setitem__('k', ConstantMechanicalLaw("k", ReducedDimension('k'), self.k))
        self.__setitem__('emissivity', ConstantMechanicalLaw("emissivity", ReducedDimension(), self.emissivity))
        self.__setitem__('E', ConstantMechanicalLaw("E", ReducedDimension('elasticModulus'), self.E))
        self.__setitem__('nu', ConstantMechanicalLaw("nu", ReducedDimension(), self.nu))
        self.__setitem__('G', ConstantMechanicalLaw("G", ReducedDimension('G'), self.G))
        self.__setitem__('alpha', ConstantMechanicalLaw("alpha", ReducedDimension(), self.alpha))
        self.__setitem__('Tref', ConstantMechanicalLaw("Tref", ReducedDimension('T'), self.Tref))
        if (self.alphaFuel is not None and self.TFuelRef is not None):
            self.__setitem__('alphaFuel', ConstantMechanicalLaw("alphaFuel", ReducedDimension(), self.alphaFuel))
            self.__setitem__('TFuelRef', ConstantMechanicalLaw("TFuelRef", ReducedDimension('T'), self.TFuelRef))
        if (self.alphaCR is not None and self.TCRRef is not None):
            self.__setitem__('alphaCR', ConstantMechanicalLaw("alphaCR", ReducedDimension(), self.alphaCR))
            self.__setitem__('TCRRef', ConstantMechanicalLaw("TCRRef", ReducedDimension('T'), self.TCRRef))
        self.__setitem__('rheologyModel', self.rheologyModel)
        return super().__repr__(depth)


    @property
    def rho(self):
        return self._rho

    @rho.setter
    def rho(self, rho) -> None:
        check_type("rho", rho, (float, int))
        self._rho = rho


    @property
    def Cp(self):
        return self._Cp

    @Cp.setter
    def Cp(self, Cp) -> None:
        check_type("Cp", Cp, (float, int))
        self._Cp = Cp


    @property
    def k(self):
        return self._k

    @k.setter
    def k(self, k) -> None:
        check_type("k", k, (float, int))
        self._k = k


    @property
    def emissivity(self):
        return self._emissivity

    @emissivity.setter
    def emissivity(self, emissivity) -> None:
        check_type("emissivity", emissivity, (float, int))
        self._emissivity = emissivity


    @property
    def E(self):
        return self._E

    @E.setter
    def E(self, E) -> None:
        check_type("E", E, (float, int))
        self._E = E


    @property
    def nu(self):
        return self._nu

    @nu.setter
    def nu(self, nu) -> None:
        check_type("nu", nu, (float, int))
        self._nu = nu


    @property
    def G(self):
        return self._G

    @G.setter
    def G(self, G) -> None:
        check_type("G", G, (float, int))
        self._G = G


    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, alpha) -> None:
        check_type("alpha", alpha, (float, int))
        self._alpha = alpha


    @property
    def Tref(self):
        return self._Tref

    @Tref.setter
    def Tref(self, Tref) -> None:
        check_type("Tref", Tref, (float, int))
        self._Tref = Tref


    @property
    def alphaFuel(self):
        return self._alphaFuel

    @alphaFuel.setter
    def alphaFuel(self, alphaFuel) -> None:
        check_type("alphaFuel", alphaFuel, (float, int), none_ok=True)
        self._alphaFuel = alphaFuel


    @property
    def TFuelRef(self):
        return self._TFuelRef

    @TFuelRef.setter
    def TFuelRef(self, TFuelRef) -> None:
        check_type("TFuelRef", TFuelRef, (float, int), none_ok=True)
        self._TFuelRef = TFuelRef


    @property
    def alphaCR(self):
        return self._alphaCR

    @alphaCR.setter
    def alphaCR(self, alphaCR) -> None:
        check_type("alphaCR", alphaCR, (float, int), none_ok=True)
        self._alphaCR = alphaCR


    @property
    def TCRRef(self):
        return self._TCRRef

    @TCRRef.setter
    def TCRRef(self, TCRRef) -> None:
        check_type("TCRRef", TCRRef, (float, int), none_ok=True)
        self._TCRRef = TCRRef


    @property
    def rheologyModel(self):
        return self._rheologyModel

    @rheologyModel.setter
    def rheologyModel(self, rheologyModel) -> None:
        check_type("rheologyModel", rheologyModel, str)
        self._rheologyModel = rheologyModel
