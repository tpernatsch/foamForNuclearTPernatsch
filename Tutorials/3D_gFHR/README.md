# 3D gFHR
Tutorial mainly prepared based on work of Yves Robert and Ludovic Jantzen
at UC Berkeley. For further details see:
Rovert, Y., et al., 2023. "IMPACT OF THERMAL COUPLING ON A PEBBLE BED 
REACTOR EQUILIBRIUM FROM HYPER-FIDELITY DEPLETION", Proceedings of the 
2023 30th International Conference on Nuclear Engineering, ICONE30
May 21-26, 2023, Kyoto, Japan. 


## Description

This is a simple 3D model of the Kairos gFHR. It show how to use the 
*nuclearSteadyStatePebble* and *lumpedNuclearStructure* power models
to describe pebbles and triso.  The *nuclearSteadyStatePebble* model is a 
steady-state model purpose made for pebble bed reactors. It allows including
several details of the pebble/triso geometry and composition in an intuitive way.
The *lumpedNuclearStructure* model is a general model that allows the user to
set-up and electric equivalent made of resistances and capacitances disposed in
series. It can be used for transient calculations but it requires the users to 
translate their problem into a electrical equvalent. For 
this tutorial, such derivation is provided in the *derivations* document, with
calculation that are done in the *lumped_structure* pythons script.

The Allrun file will simply run two different steady-states with the two models.
One can check that the two models are equivalent for instance by plotting the max
triso temperatures in paraview. 


## How to run

Use `./Allclean` to remove previous results, if any. Use `./Allrun` to
automatically run the 2 steady-states.


