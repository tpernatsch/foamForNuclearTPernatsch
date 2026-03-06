# <b>Thermal Expansion</b>

Thermal expansion defines the change in material dimensions as a function of
temperature and is a key thermo-mechanical material property in OFFBEAT.

The thermal expansion model is selected with the
<code><b>thermalExpansionModel</b></code> keyword inside the material
subdictionary of the <code><b>solverDict</b></code> (located in the
<code><b>constant</b></code> folder).

The selected model provides the thermal strain contribution induced by
temperature changes.

## <b>Classes</b>

Currently, OFFBEAT supports the following thermal expansion models:

- [constant](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionConstant.md)  
  assigns a constant thermal expansion coefficient.

- [UO2Relap](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionRelapUO2.md)  
  RELAP-based thermal expansion correlation for UO₂.

- [UPuO2Matpro](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionMatproUPuO2.md)  
  MATPRO-based correlation for (U,Pu)O₂ fuel.

- [UPuO2Martin](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionMartinUPuO2.md)  
  Martin correlation for (U,Pu)O₂ fuel.

- [UPuO2Lemehov](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionLemehovUPuO2.md)  
  Lemehov correlation for (U,Pu)O₂ fuel.

- [MaUPuO2Kato](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionMAMOX.md)  
  Kato correlation for MA-bearing MOX fuel.

- [ZircaloyMatpro](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionMatproZircaloy.md)  
  MATPRO-based correlation for Zircaloy.

- [Steel1515TiGehr](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionGehr1515Ti.md)  
  Gehr correlation for 15-15Ti steel.

- [HastelloyNSwindeman](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionSwindemanHastelloyN.md)  
  Swindeman correlation for Hastelloy N.

- [Molybdenum](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionMolybdenum.md)  
  Thermal expansion correlation for molybdenum.

- [BufferParfume](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionParfumeBuffer.md)  
  PARFUME-based correlation for buffer layers.

- [PyCParfume](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionParfumePyC.md)  
  PARFUME-based correlation for PyC layers.

- [SiCParfume](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionParfumeSiC.md)  
  PARFUME-based correlation for SiC.

- [SiCSnead](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionSneadSiC.md)  
  Snead correlation for SiC.
