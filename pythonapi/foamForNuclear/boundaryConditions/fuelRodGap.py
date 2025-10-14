from foamForNuclear.checkvalue import check_type, check_value
from .boundaryCondition import Patch


_REGION_COUPLED_TYPES = {"regionCoupledOFFBEAT"}


class FuelRodGap(Patch):
    """
    Coupled boundary condition for modeling the presence of an evolving gap
    conductance between two bodies. The gap conductance is based primarily
    on the FRAPCON 4.0 manual, making it most suitable for standard fuel rods
    (fuel pellet + cladding).

    The total gap conductance $h$ is calculated as the sum of:

    - h_gas: Gas conductance, which depends on the gap width, surface \
        roughness, jump distance and the properties of the gas mixture in the gap.
    - h_rad: Radiative heat transfer, modeled based on the emissivity of \
        the surfaces and the temperature difference.
    - h_contact: Contact conductance, which depends on the surface roughness \
        and contact pressure

    Warning: This boundary condition assumes that the same boundary condition is applied to a
    pair of `regionCoupledOFFBEAT` patches, which are coupled using the AMI
    (Arbitrary Mesh Interface) mapping algorithm.

    Warning: This boundary condition requires the presence of a `gapGasModel` in the
    simulation.

    Parameters
    ----------
    patchType : str
        Specifies the type of underlying fvPatch and must be set to
        `regionCoupledOFFBEAT`.
    kappa : str
        The name of the conductivity field (default to `k`, which is the field
        name for conductivity used by the OFFBEAT material class).
    emissivity : str
        The name of the emissivity field (default to `emissivity`, which is the
        field name for emissivity used by the OFFBEAT material class).
    gapGasModel : str
        The name of the gapGasModel class (default to `gapGas`, which is the
        name for the OFFBEAT gap gas model class).
    alpha
        Since this class is a gap-conductance model itself, it is not
        necessary to specify the gap conductance as it is calculated by the code.
        If `alpha` is specified, it will be used as the initial condition for the gap conductance.
    roughness
        Specifies the roughness field for the surface, used to calculate
        the effective gap width.
    relax : float
        Relaxation factor for the gap conductance (default to `1.0`, but values as
        low as 0.1 are often very useful to aid (and sometimes accelerate)
        convergence).
    value : float
        Specifies the initial value of the patch temperature field.
    """

    def __init__(
            self,
            value: float,
            roughness: float,
            coupled: bool=True,
            patchType: str="regionCoupledOFFBEAT",
            kappa: str="k",
            emissivity: str="emissivity",
            gapGasModel: str="gapGas",
            relax: float=1
        ):
        super().__init__(type="fuelRodGap", value=value)

        self.patchType = patchType
        self.kappa = kappa
        self.emissivity = emissivity
        self.gapGasModel = gapGasModel
        self.coupled = coupled
        self.roughness = roughness
        self.relax = relax


    def __repr__(self, depth = 0):
        self.__setitem__('patchType', self.patchType)
        self.__setitem__('kappa', self.kappa)
        self.__setitem__('emissivity', self.emissivity)
        self.__setitem__('gapGasModel', self.gapGasModel)
        self.__setitem__('coupled', self.coupled)
        self.__setitem__('relax', self.relax)
        self.__setitem__('roughness', f"uniform {self.roughness}")

        return super().__repr__(depth)


    @property
    def patchType(self):
        return self._patchType

    @patchType.setter
    def patchType(self, patchType) -> None:
        check_type("patchType", patchType, str)
        check_value("patchType", patchType, _REGION_COUPLED_TYPES)
        self._patchType = patchType


    @property
    def kappa(self):
        return self._kappa

    @kappa.setter
    def kappa(self, kappa) -> None:
        check_type("kappa", kappa, str)
        self._kappa = kappa


    @property
    def emissivity(self):
        return self._emissivity

    @emissivity.setter
    def emissivity(self, emissivity) -> None:
        check_type("emissivity", emissivity, str)
        self._emissivity = emissivity


    @property
    def gapGasModel(self):
        return self._gapGasModel

    @gapGasModel.setter
    def gapGasModel(self, gapGasModel) -> None:
        check_type("gapGasModel", gapGasModel, str)
        self._gapGasModel = gapGasModel


    @property
    def coupled(self):
        return self._coupled

    @coupled.setter
    def coupled(self, coupled) -> None:
        check_type("coupled", coupled, bool)
        self._coupled = coupled


    @property
    def relax(self):
        return self._relax

    @relax.setter
    def relax(self, relax) -> None:
        check_type("relax", relax, (int, float))
        self._relax = relax


    @property
    def roughness(self):
        return self._roughness

    @roughness.setter
    def roughness(self, roughness) -> None:
        check_type("roughness", roughness, (int, float))
        self._roughness = roughness
