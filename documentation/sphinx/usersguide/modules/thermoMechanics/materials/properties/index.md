# <b>Thermo-mechanical Properties</b>

This section documents the **thermo-mechanical material properties** available in
OFFBEAT. These properties are used by the thermal and mechanics solvers and are
defined **per material** inside the `materials` dictionary of `solverDict`.

Thermo-mechanical properties are intrinsic properties (or correlations) such as:
thermal conductivity, heat capacity, density, emissivity, thermal expansion, and
basic elastic constants.

Each property is selected through a dedicated keyword in the material definition
(e.g. `conductivityModel`, `heatCapacityModel`, …). The selected model can depend
on temperature, burnup, porosity, composition, or other state variables,
depending on the specific correlation.

## <b>Contents</b>

- [Conductivity](conductivity/index.md)
- [Heat capacity](heatCapacity/index.md)
- [Density](density/index.md)
- [Emissivity](emissivity/index.md)
- [Thermal expansion](thermalExpansion/index.md)
- [Young modulus](YoungModulus/index.md)
- [Poisson ratio](PoissonRatio/index.md)
