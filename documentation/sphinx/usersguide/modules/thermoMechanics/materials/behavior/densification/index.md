# <b>Densification</b>

Densification models describe the reduction of material volume due to pore
closure, typically occurring at the beginning of irradiation.

The densification model is selected with the
<code><b>densificationModel</b></code> keyword inside the material subdictionary
of the <code><b>solverDict</b></code>.

## <b>Classes</b>

Currently, OFFBEAT supports the following densification models:

- [none](../../../../classes/materials/materialModel/behavioralModels/densification/densificationModel.md)  
  neglects densification effects.

- [empirical](../../../../classes/materials/materialModel/behavioralModels/densification/densificationEmpirical.md)  
  empirical densification model.

- [UO2Frapcon](../../../../classes/materials/materialModel/behavioralModels/densification/densificationFrapcon.md)  
  FRAPCON-based densification model for UO₂ fuel.
