from foamForNuclear.checkvalue import check_type, check_value
from .boundaryCondition import Patch


_REGION_COUPLED_TYPES = {"regionCoupledOFFBEAT"}


class FuelRodGap(Patch):
    """
    Coupled boundary condition for modeling an **evolving thermal gap conductance**
    between two solid bodies, primarily intended for standard nuclear fuel rods
    (fuel pellet + cladding).

    The gap conductance formulation is based mainly on the FRAPCON-4.0 methodology.
    The total gap conductance :math:`h` is computed internally as the sum of:

    - **h_gas**: Gas conductance, depending on gap width, surface roughness,
      jump distance, and gap gas properties.
    - **h_rad**: Radiative heat transfer, based on surface emissivities and
      temperature difference.
    - **h_contact**: Solid-to-solid contact conductance, depending on surface
      roughness and contact pressure.

    This boundary condition must be applied symmetrically on a *pair* of
    `regionCoupledOFFBEAT` patches coupled via an AMI
    (Arbitrary Mesh Interface).

    A `gapGasModel` must be present in the simulation for this boundary
    condition to operate correctly.

    Options
    -------
    value : scalar
        Initial value of the patch temperature field.
        (required: True)

    patchType : word
        Type of the underlying fvPatch. Must be set to
        `regionCoupledOFFBEAT`.
        (default: regionCoupledOFFBEAT; required: False)

    roughness : scalar
        Surface roughness used to compute the effective gap width.
        (required: True)

    coupled : bool
        Enable coupling between the two opposing patches.
        (default: True; required: False)

    kappa : word
        Name of the thermal conductivity field.
        (default: k; required: False)

    emissivity : word
        Name of the emissivity field.
        (default: emissivity; required: False)

    gapGasModel : word
        Name of the gap gas model class used to evaluate gas properties.
        (default: gapGas; required: False)

    relax : scalar
        Relaxation factor applied to the gap conductance update.
        Values below 1 (e.g. 0.1–0.5) are often useful to improve convergence.
        (default: 1; required: False)
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
