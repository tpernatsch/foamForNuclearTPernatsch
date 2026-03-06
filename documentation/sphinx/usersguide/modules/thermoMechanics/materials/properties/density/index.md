# <b>Density</b>

Density defines the mass per unit volume of a material and is a fundamental
thermo-mechanical material property in OFFBEAT.

The density model is selected with the
<code><b>densityModel</b></code> keyword inside the material subdictionary of the
<code><b>solverDict</b></code> (located in the <code><b>constant</b></code>
folder).

Depending on the selected model, density may be treated as constant or defined
through empirical correlations.

## <b>Classes</b>

Currently, OFFBEAT supports the following density models:

- [constant](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/density/densityConstant.md)  
  assigns a constant density value.

- [UO2Constant](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/density/constantDensityUO2.md)  
  constant density for UO₂.

- [UPuO2Constant](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/density/constantDensityUPuO2.md)  
  constant density for (U,Pu)O₂ (MOX).

- [ZircaloyIaea](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/density/densityIaeaZircaloy.md)  
  IAEA correlation for Zircaloy density.

- [Steel1515TiSchumann](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/density/densitySchumann1515Ti.md)  
  Schumann correlation for 15-15Ti steel.

- [MolybdenumConstant](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/density/constantDensityMolybdenum.md)  
  constant density for molybdenum.
