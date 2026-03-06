# <b>Constitutive Laws</b>

Constitutive laws define the mechanical response of a material, i.e. how stress
is computed from the strain history and internal variables (elasticity,
plasticity, creep, and their combinations).

In OFFBEAT, the constitutive law is defined **per material** inside the
`materials` dictionary of the `solverDict` (located in the `constant` folder).

Options that are related to the general rheological framework (e.g. plane stress/strain 
assumptions and optional thermal stress settings) are set in the `rheologyOptions`
dictionary in `solverDict`. These options are documented in [rheology](../../solverDict/rheology.md).

## <b>Contents</b>

- [Laws](laws/index.md)
- [Yield Stress Models](yield_stress/index.md)
- [Creep Models](creep/index.md)
