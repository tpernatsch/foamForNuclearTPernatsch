# <b>Failure Models</b>

Failure models define criteria for material failure, loss of integrity, or
unphysical conditions during a simulation.

The failure model is selected with the
<code><b>failureModel</b></code> keyword inside the material subdictionary
of the <code><b>solverDict</b></code>.

## <b>Classes</b>

Currently, OFFBEAT supports the following failure models:

- [none](../../../../classes/materials/materialModel/behavioralModels/failureModels/failureModel.md)  
  no failure criterion applied.

- [combined](../../../../classes/materials/materialModel/behavioralModels/failureModels/combinedFailureCriteria.md)  
  combined failure criteria.

- [UO2Matpro](../../../../classes/materials/materialModel/behavioralModels/failureModels/UO2meltingMatpro.md)  
  MATPRO-based melting criterion for UO₂ fuel.

- [UPuO2Magni2020](../../../../classes/materials/materialModel/behavioralModels/failureModels/UPuO2meltingMagni2020.md)  
  Magni (2020) melting criterion for MOX fuel.

- [ZircaloyOverstrainBison](../../../../classes/materials/materialModel/behavioralModels/failureModels/ZircaloyOverstrainBison.md)  
  BISON-based overstrain failure criterion for Zircaloy.

- [ZircaloyOverstressBison](../../../../classes/materials/materialModel/behavioralModels/failureModels/ZircaloyOverstressBison.md)  
  BISON-based overstress failure criterion for Zircaloy.

- [ZircaloyPlasticInstabilityBison](../../../../classes/materials/materialModel/behavioralModels/failureModels/ZircaloyPlasticInstabilityBison.md)  
  Plastic instability criterion for Zircaloy.

- [ZircaloyRiaJernkvistModified](../../../../classes/materials/materialModel/behavioralModels/failureModels/ZircaloyRIAJernkvistModified.md)  
  Modified Jernkvist RIA failure criterion for Zircaloy.

- [ZircaloySed](../../../../classes/materials/materialModel/behavioralModels/failureModels/Zircaloy_CSED_RIA_EPRI.md)  
  CSED-RIA failure criterion for Zircaloy.
