# 2D ESFR UTOP: Steady-State Fuel Behaviour

The scope of this exercise is to simulate a Unprotected Transient OverPower and study the effects on an irradiated fuel. To achieve this, the first step is to perform the irradiation study of the fuel, starting from fresh conditions to about 600 EFPD.
Because of the complexity of the fuel-behaviour solver, this exercise is very short and is not meant to explain all the nuances of the solver, but rather to illustrate its logic and its main dictionary, i.e. constant/solverDict


---

# Tasks

1. **Complete the constant/fuelBehaviourRegion/solverDict dictionary.** In solverDict:
    - Select the correct solvers for the thermal and mechanical solutions
    - in "materials", fill the material properties to represent MOX fuel
2. **Run case.** Run the simulation with `GeN-Foam | tee log.GeN-Foam`
