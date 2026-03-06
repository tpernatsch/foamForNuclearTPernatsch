# <b>Poisson Ratio</b>

The Poisson ratio defines the ratio between transverse and axial strain under
uniaxial loading and is a basic elastic material property in OFFBEAT.

The Poisson ratio model is selected with the
<code><b>PoissonRatioModel</b></code> keyword inside the material subdictionary
of the <code><b>solverDict</b></code> (located in the <code><b>constant</b></code>
folder).

Poisson ratio models are typically constant or weakly temperature-dependent.

## <b>Classes</b>

Currently, OFFBEAT supports the following Poisson ratio models:

- [constant](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/PoissonRatioConstant.md)  
  assigns a constant Poisson ratio value.

- [UO2Constant](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/constantPoissonRatioUO2.md)  
  constant Poisson ratio for UO₂.

- [UPuO2Constant](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/constantPoissonRatioUPuO2.md)  
  constant Poisson ratio for (U,Pu)O₂ fuel.

- [ZircaloyConstant](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/constantPoissonRatioZircaloy.md)  
  constant Poisson ratio for Zircaloy.

- [ZircaloyMatpro](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/PoissonRatioMatproZircaloy.md)  
  MATPRO-based correlation for Zircaloy.

- [Steel1515TiTobbe](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/PoissonRatioTobbe1515Ti.md)  
  Tobbe correlation for 15-15Ti steel.

- [MolybdenumConstant](../../../../classes/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/constantPoissonRatioMolybdenum.md)  
  constant Poisson ratio for molybdenum.
