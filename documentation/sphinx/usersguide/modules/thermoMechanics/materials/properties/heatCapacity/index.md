# <b>Heat Capacity</b>

Heat capacity defines the amount of energy required to change the temperature of
a material and is a fundamental thermo-mechanical material property in OFFBEAT.

The heat capacity model is selected with the
<code><b>heatCapacityModel</b></code> keyword inside the material subdictionary of
the <code><b>solverDict</b></code> (located in the <code><b>constant</b></code>
folder).

The selected model provides the heat capacity as a function of temperature and,
depending on the specific correlation, may also depend on composition or other
material state variables.

## <b>Classes</b>

Currently, OFFBEAT supports the following heat capacity models:

- [constant](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityConstant.md)  
  assigns a constant heat capacity value.

- [UO2Matpro](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityMatproUO2.md)  
  MATPRO-based correlation for UO₂ fuel.

- [UPuO2Matpro](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityMatproUPuO2.md)  
  MATPRO-based correlation for (U,Pu)O₂ (MOX) fuel.

- [UPuO2Fink](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityFinkUPuO2.md)  
  Fink correlation for (U,Pu)O₂ fuel.

- [ZircaloyMatpro](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityMatproZircaloy.md)  
  MATPRO-based correlation for Zircaloy.

- [ZircaloyIaea](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityIaeaZircaloy.md)  
  IAEA correlation for Zircaloy.

- [Steel1515TiBanerjee](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityBanerjee1515Ti.md)  
  Banerjee correlation for 15-15Ti steel.

- [SiCSnead](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacitySneadSiC.md)  
  Snead correlation for SiC.

- [Molybdenum](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityMolybdenum.md)  
  Heat capacity correlation for molybdenum.
