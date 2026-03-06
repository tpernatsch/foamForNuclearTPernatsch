# **Initial and Boundary Conditions**

Setting up **initial and boundary conditions** is crucial to ensure simulations are physically meaningful and achieve proper convergence. In OFFBEAT, as in all OpenFOAM applications, these conditions are defined in the field dictionary files located in the `0/` folder (or the starting time folder). Each required field must have a corresponding file containing the initial and boundary conditions.

## **Required Fields**

Most OFFBEAT simulations require initial and boundary conditions for the following fields:

- **Temperature field (`T`)**
- **Displacement field**, which differs based on the solver:
    - `D` for total-field displacement solvers (e.g., `smallStrain`, `largeStrainTotLag`).
    - `DD` for incremental solvers (e.g., `smallStrainIncrementalUpdated`, `largeStrainUpdLag`).

Additional fields may need to be defined depending on the activated physics. For instance:

- **Neutronics solvers** require initial and boundary conditions for `neutronFlux`.
- **Element transport solvers** may involve additional field definitions.

???+ note
    
    When restarting a simulation, the code reads all field files and their initial/boundary conditions from the latest available time step. However, when starting a new simulation, only the mandatory fields (e.g., `T` and `D`) are typically required. Defining additional fields is generally unnecessary and not recommended.

This guide focuses on the primary fields `T` and `D`. Future updates will expand coverage to other fields.

___

## **Field File Structure**

An example of a typical field file for `T` is provided below:

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     9.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      T;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 1 0 0 0];

internalField   uniform 293;

boundaryField
{
    "wedge.*|.*Front|.*Back"
    {
        type            wedge;
    }

    cladOuter
    {
        type           fixedValue;
        value          uniform 500;
    }

    ".*Top.*|.*Bottom.*"
    {
        type            zeroGradient; 
    }

    etc...
}
```

The key sections in the `T` field file are:

1. **The header:**
    Contains metadata about the field, including:
    - Field name (`object` entry)
    - Field type (e.g., `volScalarField`, `volVectorField`, `volTensorField`, `volSymmTensorField`)
    - File location (e.g., `"0"`)

2. **`dimensions`:**
    Specifies the field's dimensional units in SI format (e.g., `[mass length time temperature quantity current luminousIntensity]`). These must align with OFFBEAT libraries and solvers.

3. **`internalField`:**
    Defines the field's internal values. Most cases use `uniform` values (e.g., `uniform 293` for a scalar field like temperature, or `uniform (0 0 0)` for a vector field like displacement). Non-uniform fields can be defined using OpenFOAM tools like `topoSet` or `codeStream`. Refer to the OpenFOAM documentation for details.

4. **`boundaryField`:**
    This subdictionary defines the boundary conditions for all patches in the domain. Typically, one patchField dictionary is required per patch. However, users can simplify the setup by using wildcards (e.g., `top.*`) and the `|` symbol (e.g., `top.*|bottom.*`) to group patches with similar names into a single definition. This can reduce redundancy and streamline the configuration process of the `boundaryField` subdictionary.

___

## **Boundary Conditions**

OFFBEAT uses boundary condition definitions similar to standard OpenFOAM cases. For general instructions, consult the relevant OpenFOAM user guide for your version:

- [OpenFOAM.org user guide](https://doc.cfd.direct/openfoam/user-guide-v9/)
- [OpenFOAM.com user guide](https://www.openfoam.org/documentation/user-guide)

Boundary conditions in OFFBEAT are referred to as `fvPatchField` classes in OpenFOAM terminology. Users can choose from standard OpenFOAM boundary conditions (e.g., `fixedValue`, `zeroGradient`) or use specific conditions developed in OFFBEAT for fuel performance simulations. These specialized boundary conditions are detailed in:

- [Temperature BCs](temperature/index.md)
- [Displacement BCs](displacement/index.md)


???+ warning

    For fuel-rod simulations, it is tipically required to define the gap gas composition and initial pressure in the `gapGas` file, located in the starting time folder. See the [Gap Gas Model documentation](../solverDict/gap_gas.md) for details.