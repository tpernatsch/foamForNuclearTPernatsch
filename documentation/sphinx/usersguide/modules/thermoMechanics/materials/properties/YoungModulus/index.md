# <b>Young Modulus</b>

The Young modulus defines the elastic stiffness of a material and is a key
mechanical material property in OFFBEAT.

The Young modulus model is selected with the
<code><b>YoungModulusModel</b></code> keyword inside the material subdictionary
of the <code><b>solverDict</b></code> (located in the <code><b>constant</b></code>
folder).

The selected model provides the temperature-dependent elastic modulus used by
the mechanical solver.

## <b>Classes</b>

Currently, OFFBEAT supports the following Young modulus models:

- [constant](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusConstant.md)  
  assigns a constant Young modulus value.

- [UO2Matpro](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusMatproUO2.md)  
  MATPRO-based correlation for UO₂ fuel.

- [UPuO2Matpro](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusMatproUPuO2.md)  
  MATPRO-based correlation for (U,Pu)O₂ fuel.

- [UPuO2SckCen](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusSckCenUPuO2.md)  
  SCK•CEN correlation for (U,Pu)O₂ fuel.

- [ZircaloyMatpro](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusMatproZircaloy.md)  
  MATPRO-based correlation for Zircaloy.

- [Steel1515TiTobbe](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusTobbe1515Ti.md)  
  Tobbe correlation for 15-15Ti steel.

- [SteelD9Hofman](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusHofmanD9.md)  
  Hofman correlation for D9 steel.

- [HastelloyNWatrous](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusWatrousHastelloyN.md)  
  Watrous correlation for Hastelloy N.

- [Molybdenum](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusMolybdenum.md)  
  Young modulus correlation for molybdenum.

- [BufferParfume](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusParfumeBuffer.md)  
  PARFUME-based correlation for buffer layers.

- [PyCParfume](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusParfumePyC.md)  
  PARFUME-based correlation for PyC layers.

- [SiCParfume](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusParfumeSiC.md)  
  PARFUME-based correlation for SiC.

- [SiCSnead](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusSneadSiC.md)  
  Snead correlation for SiC.
