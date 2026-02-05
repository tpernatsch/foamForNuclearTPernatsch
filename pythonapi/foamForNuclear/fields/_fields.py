from __future__ import annotations
from typing import Any, Dict

import foamlib
import numpy as np

import foamForNuclear
from foamForNuclear.mesh.blockMesh import Face
from foamForNuclear.boundaryConditions import *
from foamForNuclear.boundaryConditions.empty import Empty
from foamForNuclear.boundaryConditions.symmetry import Symmetry
from foamForNuclear.boundaryConditions.slip import Slip
from foamForNuclear.boundaryConditions.zeroGradient import ZeroGradient
from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import *
from foamForNuclear.openfoamFile import OpenFOAMFile
from foamForNuclear.mesh import Mesh, BlockMesh

_FIELD_TYPES = {"volScalarField", "volVectorField", "volTensorField"}
_GAS_PRESSURE_TYPES = {"fromModel", "fixed", "fromList"}


class Dimension:
    """
    Dimension class. Dimension objects can be multiplied or divided to create
    compound units. Default unit is no unit.

    Example
    -------
    Example of use ::

        dim1 = ffn.Dimension(default='flux')
        dim2 = ffn.Dimension(length=3)
        dim3 = dim1 * dim2
        dim4 = dim3 / dim2

    Parameters
    ----------
    default : str
        Default unit name:

        - Neutron flux : `neutronFlux` or `flux`
        - Temperature : `T` or `temperature`
        - Velocity :  `U` or `velocity`
        - Displacement: `D` or `displacement`
        - Pressure : `p` or `pressure`
        - Surface : `S`, `A`, `surface`, or `area`
        - Volume : `V` or `volume`
        - Force : `F` or `force`
        - Energy : `E` or `energy`
        - Mass density : `rho` or `massDensity`
        - Specific heat capacity : `Cp` or `specificHeatCapacity`
        - Thermal conductivity : `k` or `thermalConductivity`
        - Turbulent kinetic energy : `turbulentKineticEnergy`
        - Turbulent kinetic energy dissipation rate : `epsilon` or `turbulentKineticEnergyDissipationRate`
        - Specific dissipation rate : `omega` or `specificDissipationRate`
        - Turbulent (eddy) kinematic viscosity : `nut` or `turbulentKinematicViscosity`
        - Turbulent (eddy) dynamic viscosity : `mut` or `turbulentDynamicViscosity`
        - Thermal eddy diffusivity : `alphat` or `thermalEddyDiffusivity`

    mass : float
        Mass in kilogram
    length : float
        Length in metre
    time : float
        Time in seconds
    temperature : float
        Temperature in Kelvin
    moles : float
        Moles in mole
    current : float
        Current in Ampere
    luminousyIntensity : float
        Luminous intensity in Candela

    Attributes
    ----------
    mass : float
        Mass in kilogram
    length : float
        Length in metre
    time : float
        Time in seconds
    temperature : float
        Temperature in Kelvin
    moles : float
        Moles in mole
    current : float
        Current in Ampere
    luminousyIntensity : float
        Luminous intensity in Candela
    """
    def __init__(
            self,
            default: str='',
            mass=0,
            length=0,
            time=0,
            temperature=0,
            moles=0,
            current=0,
            luminousyIntensity=0,
        ):
        self.mass = mass
        self.length = length
        self.time = time
        self.temperature = temperature
        self.moles = moles
        self.current = current
        self.luminousyIntensity = luminousyIntensity

        if (default in ["flux", "neutronFlux"]):
            self.set_to_flux()
        elif (default in ["T", "temperature"]):
            self.set_to_temperature()
        elif (default in ["U", "velocity"]):
            self.set_to_velocity()
        elif (default in ["D", "displacement"]):
            self.set_to_displacement()
        elif (default in ["p", "pressure", "elasticModulus", "G"]):
            self.set_to_pressure()
        elif (default in ["S", "A", "surface", "area"]):
            self.set_to_surface()
        elif (default in ["V", "volume"]):
            self.set_to_volume()
        elif (default in ["F", "force"]):
            self.set_to_force()
        elif (default in ["E", "energy"]):
            self.set_to_energy()
        elif (default in ["rho", "massDensity"]):
            self.set_to_mass_density()
        elif (default in ["Cp", "specificHeatCapacity"]):
            self.set_to_specific_heat_capacity()
        elif (default in ["k", "thermalConductivity"]):
            self.set_to_thermal_conductivity()
        elif (default in ["turbulentKineticEnergy"]):
            self.set_to_turbulent_kinetic_energy()
        elif (default in ["epsilon", "turbulentKineticEnergyDissipationRate"]):
            self.set_to_turbulent_kinetic_energy_dissipation_rate()
        elif (default in ["omega", "specificDissipationRate"]):
            self.set_to_specific_dissipation_rate()
        elif (default in ["nut", "turbulentKinematicViscosity"]):
            self.set_to_turbulent_kinematic_viscosity()
        elif (default in ["mut", "turbulentDynamicViscosity"]):
            self.set_to_turbulent_dynamic_viscosity()
        elif (default in ["alphat", "thermalEddyDiffusivity"]):
            self.set_to_thermal_eddy_diffusivity()


    def __repr__(self):
        return(f"[{self.mass} {self.length} {self.time} {self.temperature} {self.moles} {self.current} {self.luminousyIntensity}]")

    def __mul__(self, rhs):
        return(Dimension(
            mass=self.mass + rhs.mass,
            length=self.length + rhs.length,
            time=self.time + rhs.time,
            temperature=self.temperature + rhs.temperature,
            moles=self.moles + rhs.moles,
            current=self.current + rhs.current,
            luminousyIntensity=self.luminousyIntensity + rhs.luminousyIntensity,
        ))

    def __truediv__(self, rhs):
        return(Dimension(
            mass=self.mass - rhs.mass,
            length=self.length - rhs.length,
            time=self.time - rhs.time,
            temperature=self.temperature - rhs.temperature,
            moles=self.moles - rhs.moles,
            current=self.current - rhs.current,
            luminousyIntensity=self.luminousyIntensity - rhs.luminousyIntensity,
        ))

    def __eq__(self, rhs):
        return(
            self.mass == rhs.mass and
            self.length == rhs.length and
            self.time == rhs.time and
            self.temperature == rhs.temperature and
            self.moles == rhs.moles and
            self.current == rhs.current and
            self.luminousyIntensity == rhs.luminousyIntensity
        )

    def __ne__(self, rhs):
        return(not self.__eq__(rhs))


    def reset(self):
        self.mass = 0
        self.length = 0
        self.time = 0
        self.temperature = 0
        self.moles = 0
        self.current = 0
        self.luminousyIntensity = 0

    def set_to_flux(self):
        self.reset()
        self.length = -2
        self.time = -1

    def set_to_temperature(self):
        self.reset()
        self.temperature = 1

    def set_to_velocity(self):
        self.reset()
        self.length = 1
        self.time = -1

    def set_to_displacement(self):
        self.reset()
        self.length = 1

    def set_to_pressure(self):
        self.reset()
        self.mass = 1
        self.length = -1
        self.time = -2

    def set_to_surface(self):
        self.reset()
        self.length = 2

    def set_to_volume(self):
        self.reset()
        self.length = 3

    def set_to_force(self):
        self.reset()
        self.mass = 1
        self.length = 1
        self.time = -2

    def set_to_energy(self):
        self.reset()
        self.mass = 1
        self.length = 2
        self.time = -2

    def set_to_mass_density(self):
        self.reset()
        self.mass = 1
        self.length = -3

    def set_to_specific_heat_capacity(self):
        self.reset()
        self.length = 2
        self.time = -2
        self.temperature = -1

    def set_to_thermal_conductivity(self):
        self.reset()
        self.mass = 1
        self.length = 1
        self.time = 3
        self.temperature = -1

    def set_to_turbulent_kinetic_energy(self):
        self.reset()
        self.length = 2
        self.time = -2

    def set_to_turbulent_kinetic_energy_dissipation_rate(self):
        self.reset()
        self.length = 2
        self.time = -3

    def set_to_specific_dissipation_rate(self):
        self.reset()
        self.time = -1

    def set_to_turbulent_kinematic_viscosity(self):
        self.reset()
        self.length = 2
        self.time = -1

    def set_to_turbulent_dynamic_viscosity(self):
        self.reset()
        self.mass = 1
        self.length = -1
        self.time = -1

    def set_to_thermal_eddy_diffusivity(self):
        self.reset()
        self.mass = 1
        self.length = -1
        self.time = -1


class ReducedDimension(Dimension):
    """
    Reduced dimension object used for material declaration in thermo-mechanics
    solver.
    """
    def __init__(
            self,
            default: str='',
            mass: int=0,
            length: int=0,
            time: int=0,
            temperature: int=0,
            moles: int=0
        ):
        super().__init__(default, mass, length, time, temperature, moles)

    def __repr__(self):
        return(f"[{self.mass} {self.length} {self.time} {self.temperature} {self.moles}]")


class Field(OpenFOAMFile):
    """
    Field object collecting the internal value per cell and the boundary
    conditions.

    Parameters
    ----------
    name : str
        Name of the field
    dimensions : Dimension
        Dimension of the field (default to no unit)
    internalField : int | float | Vector
        Value of the internal field
    boundaryField : dict | OpenFOAMDict
        Dictionary grouping all the boundary conditions attached to the field
    region : str
        Name of the region (default to `""`).
    """
    def __init__(
            self,
            name: str,
            dimensions: Dimension=Dimension(),
            internalField: int | float | list[float] | np.ndarray[float] | Vector=None,
            boundaryField: dict | OpenFOAMDict=None,
            region: str="",
            default_bc: Patch = ZeroGradient()
        ):
        super().__init__(name, region=region)

        self.fieldType = None
        self.dimensions = dimensions
        self.internalField = internalField
        self.boundaryField = boundaryField
        self.default_bc = default_bc
        self._explicit_bc_keys: set[str] = set()

    @property
    def dimensions(self):
        return self._dimensions

    @dimensions.setter
    def dimensions(self, dimensions) -> None:
        check_type("dimensions", dimensions, Dimension)
        self._dimensions = dimensions

    @property
    def internalField(self):
        return self._internalField

    @internalField.setter
    def internalField(self, internalField) -> None:
        check_type("internalField", internalField, (float, int, list, np.ndarray, Vector), none_ok=True)
    
        if (isinstance(internalField, (list, np.ndarray))
                and len(internalField) == 3):
            self._internalField = Vector(internalField[0], internalField[1], internalField[2])
        else:
            self._internalField = internalField
        if (internalField is not None):
            self.setFieldType()

    @property
    def boundaryField(self):
        return self._boundaryField

    @boundaryField.setter
    def boundaryField(self, boundaryField) -> None:
        check_type("boundaryField", boundaryField, (OpenFOAMDict, dict), none_ok=True)
        if (boundaryField is not None):
            self._boundaryField = OpenFOAMDict(boundaryField)
        else:
            self._boundaryField = OpenFOAMDict({}) # Do not put it as default in the __init__

    @property
    def default_bc(self):
        return self._default_bc

    @default_bc.setter
    def default_bc(self, default_bc) -> None:
        check_type("default_bc", default_bc, Patch)
        self._default_bc = default_bc

    def getPatch(self, dict_):
        if (dict_['type'] == "calculated"):
            return(Calculated())
        elif (dict_['type'] == "empty"):
            return(Empty())
        elif (dict_['type'] == "fixedValue"):
            return(FixedValue(value=dict_['value']))
        elif (dict_['type'] == "slip"):
            return(Slip())
        elif (dict_['type'] == "wedge"):
            return(Wedge())
        elif (dict_['type'] == "zeroGradient"):
            return(ZeroGradient())
        # To be continued
        return(Empty())


    def setFieldType(self):
        if (isinstance(self.internalField, (int, float))):
            self.fieldType = 'volScalarField'
        elif (isinstance(self.internalField, (list, Vector, np.ndarray))):
            self.fieldType = 'volVectorField'
        elif (isinstance(self.internalField, (list, Tensor, np.ndarray))):
            self.fieldType = 'volTensorField'
        else:
            msg = f"Cannot find a proper field format for type {type(self.internalField).__name__}"
            raise ValueError(msg)

    def prepare_boundary_field(
        self,
        mesh: BlockMesh,
        *,
        overwrite: bool = False,
        auto_defaults: bool = True
    ) -> list[str]:
        """
        Takes `mesh` and ensure `self.boundaryField` has one entry per patch.

        Strategy:
        - Track concrete patch names explicitly set by the user in `self._set_patch_names`.
        - Keep regex/wildcard entries collapsed in `boundaryField` and treat them as covering
            the matching real patches (so we don't auto-fill those patches).
        - Auto-fill defaults ONLY for patches that are not explicitly set AND not covered by any regex key.

        Returns the list of patch names it touched (concrete patch names only).
        """
        import re
        from typing import Any

        # ---- get patches ----
        patches = [
            face for face in mesh.faces
            if face.name != "pipeWall" and not getattr(face, "toBeMerged", False)
        ]

        # ---- ensure storage ----
        if not hasattr(self, "boundaryField") or self.boundaryField is None:
            self.boundaryField = {}

        # ---- ensure explicit-set tracking ----
        # This is the "set keys" you want.
        if not hasattr(self, "_set_patch_names") or self._set_patch_names is None:
            self._set_patch_names: set[str] = set()

        # ---- helper: pick a default BC for a given patch ----
        def _default_for(patch) -> Any:
            if not auto_defaults:
                return self.default_bc
            try:
                kind = patch.boundaryType
            except Exception:
                return self.default_bc
            if kind == "wedge":
                return Wedge()
            if kind == "empty":
                return Empty()
            if kind == "symmetry":
                return Symmetry()
            if kind == "cyclic":
                return Cyclic()
            if kind == "cyclicAMI":
                return CyclicAMI()
            return self.default_bc

        # ---- detect regex keys (collapsed entries) ----
        # You hinted: "they always start with \"".
        # We'll support that explicitly and also fall back to a light heuristic.
        def _is_regex_key(k: str) -> bool:
            k_strip = k.strip()
            if k_strip.startswith('"') and k_strip.endswith('"') and len(k_strip) >= 2:
                return True
            # fallback heuristic for unquoted regex keys
            return any(ch in k for ch in ["|", "*", "(", ")", "[", "]", "^", "$", "+", "?", "."])

        def _normalize_regex_key(k: str) -> str:
            k = k.strip()
            # remove surrounding quotes if present
            if k.startswith('"') and k.endswith('"') and len(k) >= 2:
                k = k[1:-1]
            # remove spaces (so "a | b" works)
            k = "".join(k.split())
            return k

        # collect regex keys currently present in boundaryField
        existing_keys = list(self.boundaryField.keys())
        regex_keys = [k for k in existing_keys if _is_regex_key(k)]

        # precompile regex patterns once
        compiled_regex: list[re.Pattern] = []
        for k in regex_keys:
            pat = _normalize_regex_key(k)
            try:
                compiled_regex.append(re.compile(pat))
            except re.error:
                # If a user stored something that looks like regex but isn't valid python re,
                # we just ignore it for coverage (but we keep the entry in boundaryField).
                # You can raise instead if you want strict behavior.
                continue

        def _covered_by_regex(patch_name: str) -> bool:
            for rx in compiled_regex:
                if rx.fullmatch(patch_name) or rx.match(patch_name):
                    return True
            return False

        # ---- seed entries ----
        touched: list[str] = []

        for p in patches:
            pname = p.name

            if not overwrite:
                # If user explicitly set this concrete patch name, keep it.
                # (Even if it equals the default BC.)
                if pname in self._set_patch_names:
                    continue

                # If this patch is covered by any regex BC entry, don't auto-fill it.
                if _covered_by_regex(pname):
                    continue

                # If the patch exists literally in boundaryField, keep it as well.
                # (This can happen if you loaded from file or user set via dict directly.)
                if pname in self.boundaryField.keys():
                    continue

            patch_bc = _default_for(p)

            if patch_bc is not None:
                # IMPORTANT: write directly to boundaryField to avoid turning the concrete name into a regex
                self.boundaryField[pname] = patch_bc
            else:
                self.boundaryField[pname] = {}

            touched.append(pname)

        return touched


    def add_boundary_condition(self, name: str | Face, boundaryCondition: Patch):
        """
        Add/Set a boundary condition to a boundary face

        Parameters
        ----------
        name : str | Face
            Name of the patch or Face to apply the boundary condition
        boundaryCondition : Patch
            Boundary condition definition
        """
        self.set_boundary_condition(name=name, boundaryCondition=boundaryCondition)

    def set_boundary_condition(self, name: str | Face | list[str] | list[Face], boundaryCondition: Patch):
        """
        Add/Set a boundary condition to a boundary face or to a list of boundary faces

        Parameters
        ----------
        name : str | Face
            Name (or list of names) of the patch(es) or Face(s) to which the boundary condition is applied
        boundaryCondition : Patch
            Boundary condition definition
        """
        if not hasattr(self, "_explicit_bc_keys") or self._explicit_bc_keys is None:
            self._explicit_bc_keys = set()

        names = name if isinstance(name, list) else [name]

        for name in names:
            check_type("name", name, (str, Face))
            check_type("boundaryCondition", boundaryCondition, Patch)
            if (isinstance(name, str)):
                name = format_to_openfoam_regex(name)
                self.boundaryField[name] = boundaryCondition
                self._explicit_bc_keys.add(name)
            elif (isinstance(name, Face)):
                self.boundaryField[name.name] = boundaryCondition
                self._explicit_bc_keys.add(name.name)
            else:
                msg = f"Boundary name must be a 'str' OR 'Face'. Found type {type(name).__name__}"
                raise ValueError(msg)


    def export_to_openfoam(self):
        with open(f"{self.folder}/{self.region}/{self.name}", 'w') as f:
            f.write(openfoamHeader)
            f.write(openfoamFileHeader(self.name, self.fieldType))

            f.write(addParameter('dimensions', self.dimensions, isAddExtraLine=True))

            if (isinstance(self.internalField, (float, int, Vector))):
                f.write(addParameter('internalField', f"uniform {self.internalField}", isAddExtraLine=True))

            f.write(f"boundaryField{self.boundaryField!r}\n")

            f.write(openfoamFooterLine)


    def import_from_openfoam(self, path: str=None):
        foamField = foamlib.FoamFieldFile(path if path is not None else self.path)

        # Internal field
        internalField = foamField.internal_field
        if (isinstance(internalField, (list, np.ndarray))):
            # Field is a list
            if (len(internalField.shape) == 1):
                if (len(internalField) == 3):
                    self.internalField = Vector(internalField[0], internalField[1], internalField[2])
                # elif (len(internalField) == 6):
                #     self.internalField = Tensor(internalField)
            elif (len(internalField.shape) == 2):
                None
        elif (isinstance(internalField, (float, int))):
            self.internalField = internalField

        # Boundary
        self.boundaryField = OpenFOAMDict({})
        for bcName, bc in foamField.boundary_field.items():
            self.set_boundary_condition(bcName, self.getPatch(bc))

        # Dimensions
        dim = foamField.dimensions
        self.dimensions = Dimension(
            mass=dim.mass,
            length=dim.length,
            time=dim.time,
            temperature=dim.temperature,
            moles=dim.moles,
            current=dim.current,
            luminousyIntensity=dim.luminous_intensity
        )


class GapGas(OpenFOAMFile):
    """
    The `FRAPCON` class in `OFFBEAT` requires the user to set various parameters
    in `gapGasOptions` (a subdictionary of the `solverDict`) and `gapGas`
    (a dictionary in the `0/uniform/` folder), as described below.

    Parameters
    ----------
    gasPressureType : str
        Defines the method for gas pressure calculation (options: `fromModel`,
        `fixed`, or `fromList`).
    gasPressure : float
        Fixed gas pressure value when `gasPressureType` is set to `fixed`.
    gapPressureList : Table
        Activated only if `gasPressureType` is set as `fromList`, this is a
        sub-dictionary used to provide a time-dependent table for gas pressure.
    type
        Parameter of `gapPressureList` subdictionary; it specifies the table
        type (`table`).
    values
        Parameter of `gapPressureList` subdictionary; it is a list of
        time-pressure pairs (e.g., `[(0.0, 1e5), (100.0, 1e6)]`).
    outOfBounds
        Parameter of `gapPressureList` subdictionary; it handles the method for
        out-of-bounds values (default: `clamp`).
    interpolationScheme
        Parameter of `gapPressureList` subdictionary; it specifies the
        interpolation method between time points (e.g., `linear`).
    region : str
        Name of the region (default `""`)

    Attributes
    ----------
    massFractions : OpenFOAMDict
        Specifies the initial mass fractions of gases (`Ar`, `He`, `Kr`, `Ne`,
        `Rn`, `Xe`). Default to pure helium.
    """

    def __init__(
            self,
            gasPressureType="fromModel",
            gasPressure=None,
            # gapPressureList=None,
            region = ""
        ):
        super().__init__("gapGas", region=region)

        self.gasPressureType = gasPressureType
        self.gasPressure = gasPressure
        # self.gapPressureList = gapPressureList
        self.massFractions = OpenFOAMDict({
            'Ar': 0,
            'He': 1,
            'Kr': 0,
            'Ne': 0,
            'Rn': 0,
            'Xe': 0,
        })

    @property
    def gasPressureType(self):
        return self._gasPressureType

    @gasPressureType.setter
    def gasPressureType(self, gasPressureType) -> None:
        check_type("gasPressureType", gasPressureType, str)
        check_value("gasPressureType", gasPressureType, _GAS_PRESSURE_TYPES)
        self._gasPressureType = gasPressureType

    @property
    def gasPressure(self):
        return self._gasPressure

    @gasPressure.setter
    def gasPressure(self, gasPressure) -> None:
        check_type("gasPressure", gasPressure, (float, int), none_ok=True)
        self._gasPressure = gasPressure

    # @property
    # def gapPressureList(self):
    #     return self._gapPressureList

    # @gapPressureList.setter
    # def gapPressureList(self, gapPressureList) -> None:
    #     check_type("gapPressureList", gapPressureList, (Table, list, np.ndarray), none_ok=True)
    #     if (isinstance(gapPressureList, Table)):
    #         self._gapPressureList = gapPressureList
    #     elif (gapPressureList is None):
    #         self._gapPressureList = None
    #     else:
    #         self._gapPressureList = Table(gapPressureList)


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = ""

        text += f"massFractions{self.massFractions!r}\n"

        text += addParameter('gasPressureType', self.gasPressureType, isAddExtraLine=True)

        if (self.gasPressure is not None):
            text += addParameter('gasPressure', self.gasPressure, isAddExtraLine=True)

        # if (self.gapPressureList is not None):
        #     text += addParameter('gapPressureList', self.gapPressureList, isAddExtraLine=True)

        return(text)


# OFFBEAT Fields
class Temperature(Field):
    def __init__(
            self,
            internalField: int | float =None,
            boundaryField: dict | OpenFOAMDict=None,
            region: str="",
            default_bc: Patch = ZeroGradient()
        ):
        super().__init__(
            name="T",
            dimensions=Dimension(default="T"),
            internalField=internalField,
            boundaryField=boundaryField,
            region=region,
            default_bc=default_bc
        )


class NeutronFlux0(Field):
    def __init__(
            self,
            internalField: int | float =None,
            boundaryField: dict | OpenFOAMDict=None,
            region: str="",
            default_bc: Patch = ZeroGradient()
        ):
        super().__init__(
            name="neutronFlux0",
            dimensions=Dimension(default="flux"),
            internalField=internalField,
            boundaryField=boundaryField,
            region=region,
            default_bc=default_bc
        )

class Burnup(Field):
    def __init__(
            self,
            internalField: int | float =None,
            boundaryField: dict | OpenFOAMDict=None,
            region: str="",
            default_bc: Patch = ZeroGradient()
        ):
        super().__init__(
            name="Bu",
            dimensions=Dimension(default=""),
            internalField=internalField,
            boundaryField=boundaryField,
            region=region,
            default_bc=default_bc
        )


class Displacement(Field):
    def __init__(
            self,
            internalField: Vector=None,
            boundaryField: dict | OpenFOAMDict=None,
            region: str="",
            default_bc: Patch = ZeroGradient()
        ):
        super().__init__(
            name="D",
            dimensions=Dimension(default="D"),
            internalField=internalField,
            boundaryField=boundaryField,
            region=region,
            default_bc=default_bc
        )


# ffn Fields
class Tliquid(Field):
    def __init__(
            self,
            internalField: int | float =None,
            boundaryField: dict | OpenFOAMDict=None,
            region: str="",
            default_bc: Patch = ZeroGradient()
        ):
        super().__init__(
            name="T.liquid",
            dimensions=Dimension(default="T"),
            internalField=internalField,
            boundaryField=boundaryField,
            region=region,
            default_bc=default_bc
        )

class Tstructure(Field):
    def __init__(
            self,
            internalField: int | float =None,
            boundaryField: dict | OpenFOAMDict=None,
            region: str="",
            default_bc: Patch = ZeroGradient()
        ):
        super().__init__(
            name="T.structure",
            dimensions=Dimension(default="T"),
            internalField=internalField,
            boundaryField=boundaryField,
            region=region,
            default_bc=default_bc
        )

class Tvapour(Field):
    def __init__(
            self,
            internalField: int | float =None,
            boundaryField: dict | OpenFOAMDict=None,
            region: str="",
            default_bc: Patch = ZeroGradient()
        ):
        super().__init__(
            name="T.vapour",
            dimensions=Dimension(default="T"),
            internalField=internalField,
            boundaryField=boundaryField,
            region=region,
            default_bc=default_bc
        )

class Uliquid(Field):
    def __init__(
            self,
            internalField: Vector  =None,
            boundaryField: dict | OpenFOAMDict=None,
            region: str="",
            default_bc: Patch = ZeroGradient()
        ):
        super().__init__(
            name="U.liquid",
            dimensions=Dimension(default="U"),
            internalField=internalField,
            boundaryField=boundaryField,
            region=region,
            default_bc=default_bc
        )

class Uvapour(Field):
    def __init__(
            self,
            internalField: Vector =None,
            boundaryField: dict | OpenFOAMDict=None,
            region: str="",
            default_bc: Patch = ZeroGradient()
        ):
        super().__init__(
            name="U.vapour",
            dimensions=Dimension(default="U"),
            internalField=internalField,
            boundaryField=boundaryField,
            region=region,
            default_bc=default_bc
        )

class alphaLiquid(Field):
    def __init__(
            self,
            internalField: int | float =None,
            boundaryField: dict | OpenFOAMDict=None,
            region: str="",
            default_bc: Patch = ZeroGradient()
        ):
        super().__init__(
            name="alpha.liquid",
            internalField=internalField,
            boundaryField=boundaryField,
            region=region,
            default_bc=default_bc
        )

class alphaVapour(Field):
    def __init__(
            self,
            internalField: int | float =None,
            boundaryField: dict | OpenFOAMDict=None,
            region: str="",
            default_bc: Patch = ZeroGradient()
        ):
        super().__init__(
            name="alpha.vapour",
            internalField=internalField,
            boundaryField=boundaryField,
            region=region,
            default_bc=default_bc
        )

class p(Field):
    def __init__(
            self,
            internalField: int | float =None,
            boundaryField: dict | OpenFOAMDict=None,
            region: str="",
            default_bc: Patch = ZeroGradient()
        ):
        super().__init__(
            name="p",
            dimensions=Dimension(default="p"),
            internalField=internalField,
            boundaryField=boundaryField,
            region=region,
            default_bc=default_bc
        )

class p_rgh(Field):
    def __init__(
            self,
            internalField: int | float =None,
            boundaryField: dict | OpenFOAMDict=None,
            region: str="",
            default_bc: Patch = ZeroGradient()
        ):
        super().__init__(
            name="p_rgh",
            dimensions=Dimension(default="p"),
            internalField=internalField,
            boundaryField=boundaryField,
            region=region,
            default_bc=default_bc
        ) 