# **Displacement Boundary Conditions**

This section summarizes the main displacement boundary conditions currently available in OFFBEAT. Many other boundary conditions (including standard `fixedValue` and `fixedGradient`) are available in OpenFOAM and are not treated here. We refer to the OpenFOAM documentation for more information.

## **Fixed Displacement (Dirichlet) Boundary Conditions**

These boundary conditions impose a displacement value (or time series) on a patch.

Available fixed value boundary conditions include:

- [fixedDisplacement](../../../classes/fvPatchFields/fixedDisplacement/fixedDisplacementFvPatchVectorField.md) 

## **Zero Shear (Or Sliding Support) Boundary Conditions**

These boundary conditions constrain displacement along certain directions while allowing motion in others, often used to simulate roller-like supports or unilateral contact.

Available options include:

- [fixedDisplacementZeroShear](../../../classes/fvPatchFields/fixedDisplacementZeroShear/fixedDisplacementZeroShearFvPatchVectorField.md)
- [unilateralContact](../../../classes/fvPatchFields/fixedDisplacementZeroShear/unilateralContactFvPatchVectorField.md)

## **Fixed Traction or Fixed Pressure Boundary Conditions**

These boundary conditions impose fixed tractions or pressures on a patch.

Available fixed traction boundary conditions include:

- [tractionDisplacement](../../../classes/fvPatchFields/tractionDisplacement/tractionDisplacementFvPatchVectorField.md)
- [coolantPressure](../../../classes/fvPatchFields/tractionDisplacement/coolantPressureFvPatchVectorField.md)
- [gapPressure](../../../classes/fvPatchFields/tractionDisplacement/gapPressureFvPatchVectorField.md)
- [plenumSpringPressure](../../../classes/fvPatchFields/tractionDisplacement/plenumSpringPressureFvPatchVectorField.md)
- [topCladRingPressure](../../../classes/fvPatchFields/tractionDisplacement/topCladRingPressureFvPatchVectorField.md)
- [gapContact](../../../classes/fvPatchFields/tractionDisplacement/gapContactFvPatchVectorField.md)
- [fixedPointLoad](../../../classes/fvPatchFields/tractionDisplacement/fixedPointLoadFvPatchVectorField.md)
- [fixedTorque](../../../classes/fvPatchFields/tractionDisplacement/fixedTorqueFvPatchVectorField.md)


## **Implicit Contact**

These boundary conditions imposes a contact constraint between two bodies in an implicit manner.

Available implicit contact boundary conditions include:

- [implicitGapContact](../../../classes/fvPatchFields/implicitContact/implicitGapContactFvPatchVectorField.md)