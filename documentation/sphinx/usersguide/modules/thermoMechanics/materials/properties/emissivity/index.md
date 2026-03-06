# <b>Emissivity</b>

Emissivity defines the efficiency of a material surface in emitting thermal
radiation and is used by radiative heat transfer models in OFFBEAT.

The emissivity model is selected with the
<code><b>emissivityModel</b></code> keyword inside the material subdictionary of
the <code><b>solverDict</b></code> (located in the <code><b>constant</b></code>
folder).

Emissivity models are typically constant or empirical correlations depending on
material type and temperature.

## <b>Classes</b>

Currently, OFFBEAT supports the following emissivity models:

- [constant](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/emissivity/emissivityConstant.md)  
  assigns a constant emissivity value.

- [UO2Relap](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/emissivity/emissivityRelapUO2.md)  
  RELAP-based emissivity correlation for UO₂ fuel.

- [ZircaloyConstant](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/emissivity/constantEmissivityZircaloy.md)  
  constant emissivity for Zircaloy.

- [MolybdenumConstant](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/emissivity/constantEmissivityMolybdenum.md)  
  constant emissivity for molybdenum.
