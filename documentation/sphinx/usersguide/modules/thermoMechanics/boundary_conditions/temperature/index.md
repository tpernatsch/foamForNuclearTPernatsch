# **Temperature Boundary Conditions**

This section lists the main thermal boundary conditions currently available in OFFBEAT. Many other boundary conditions (including standard `fixedValue` and `fixedGradient`) are available in OpenFOAM and are not treated here. We refer to the OpenFOAM documentation for more information.

## **Fixed temperature (Dirichlet) boundary conditions**

The `fixedTemperature` boundary condition base class is similar to `fixedValue` in OpenFOAM, with the added capability of incorporating resistance for the oxide layer (when used for the outer cladding). The derived classes provide additional flexibility by enabling radial or axial time-dependent temperature profiles.

Available coupled boundary conditions include:

- [fixedTemperature](../../../classes/fvPatchFields/fixedTemperature/fixedTemperatureFvPatchScalarField.md)
- [timeDependentTemperature](../../../classes/fvPatchFields/fixedTemperature/timeDependentTemperatureFvPatchScalarField.md)
- [axialProfileT](../../../classes/fvPatchFields/fixedTemperature/axialProfileFvPatchScalarField.md)
- [timeDependentAxialProfileT](../../../classes/fvPatchFields/fixedTemperature/timeDependentAxialProfileFvPatchScalarField.md)
- [xzTemperatureProfile](../../../classes/fvPatchFields/fixedTemperature/xzTemperatureProfileFvPatchScalarField.md)

## **Coupled Boundary Conditions**

These boundary conditions are designed to handle heat exchange between two bodies, whether in direct contact or separated by a small gap. They are particularly suited for scenarios involving fuel and cladding heat exchange.

Available coupled boundary conditions include:

- [temperatureCoupled](../../../classes/fvPatchFields/temperatureCoupled/temperatureCoupledFvPatchScalarField.md)
- [resistiveGap](../../../classes/fvPatchFields/temperatureCoupled/resistiveGapFvPatchScalarField.md)
- [fuelRodGap](../../../classes/fvPatchFields/temperatureCoupled/fuelRodGapFvPatchScalarField.md)
- [trisoGap](../../../classes/fvPatchFields/temperatureCoupled/trisoGapFvPatchScalarField.md)

## **Convective Boundary Conditions**

These boundary conditions simulate convective heat exchange mechanism between a solid and a fluid. They are mostly used in OFFBEAT for cladding outer surface.

Available options include:

- [convectiveHTC](../../../classes/fvPatchFields/convectiveHTC/convectiveHTCFvPatchScalarField.md) 
- [timeDependentAxialProfileHTC](../../../classes/fvPatchFields/convectiveHTC/timeDependentAxialProfileHTCfvPatchScalarField.md)

## **Radiative-Convective Sink Boundary Conditions**

These boundary conditions simulate convective heat exchange mechanism between a solid and an environment via radiation and convection.

Available options include:

- [radiativeConvectiveSink](../../../classes/fvPatchFields/radiativeConvectiveSink/radiativeConvectiveSinkFvPatchScalarField.md)


## **Coolant channel boundary conditions**

These boundary conditions simulate convective heat exchange between a solid and a coolant flowing on its surface.

Available options include:

- [coolantChannelfvPatchScalarField](../../../classes/fvPatchFields/coolantChannel/coolantChannelfvPatchScalarField.md)
- [NaCoolantChannelfvPatchScalarField](../../../classes/fvPatchFields/coolantChannel/NaCoolantChannelfvPatchScalarField.md)

