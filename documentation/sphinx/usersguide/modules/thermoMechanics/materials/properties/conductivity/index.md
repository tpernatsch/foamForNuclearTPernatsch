# <b>Thermal Conductivity</b>

Thermal conductivity defines the ability of a material to conduct heat and is a
core material property used by the thermal solver in OFFBEAT.

The thermal conductivity model is selected with the
<code><b>conductivityModel</b></code> keyword inside the material subdictionary of
the <code><b>solverDict</b></code> (located in the <code><b>constant</b></code>
folder).

The selected model provides the thermal conductivity as a function of
temperature and, depending on the specific correlation, may also depend on
other state variables such as burnup, porosity, or composition.

## <b>Classes</b>

Currently, OFFBEAT supports the following thermal conductivity models:

- [constant](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityConstant.md)  
  assigns a constant thermal conductivity value.

- [UO2Matpro](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityMatproUO2.md)  
  MATPRO-based correlation for UO₂ fuel.

- [UO2Ifa601](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityUO2IFA601.md)  
  IFA-601 correlation for UO₂ fuel.

- [UO2Nfir](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityNfirUO2.md)  
  NFIR-based correlation for UO₂ fuel.

- [MaUPuO2Magni](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityUPuO2MagniMAMOX.md)  
  Magni correlation for (U,Pu)O₂ (MOX) fuel.

- [UPuO2Brancheria](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityUPuO2Brancheria.md)  
  Brancheria correlation for (U,Pu)O₂ fuel.

- [UPuO2Kato](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityUPuO2Kato.md)  
  Kato correlation for (U,Pu)O₂ fuel.

- [UPuO2LanningBeyer](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityUPuO2LanningBeyer.md)  
  Lanning–Beyer correlation for (U,Pu)O₂ fuel.

- [UPuO2Philipponneau](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityUPuO2Philipponneau.md)  
  Philipponneau correlation for (U,Pu)O₂ fuel.

- [ZircaloyRelap](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityRelapZircaloy.md)  
  RELAP-based correlation for Zircaloy.

- [Steel1515TiTobbe](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityTobbe1515Ti.md)  
  Tobbe correlation for 15-15Ti steel.

- [HastelloyNSwindeman](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivitySwindemanHastelloyN.md)  
  Swindeman correlation for Hastelloy N.

- [Molybdenum](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityMolybdenum.md)  
  Thermal conductivity correlation for molybdenum.

- [BufferParfume](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityParfumeBuffer.md)  
  PARFUME-based correlation for buffer layers.

- [SiCParfume](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityParfumeSiC.md)  
  PARFUME-based correlation for SiC.
