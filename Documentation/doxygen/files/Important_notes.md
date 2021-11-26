# Important notes {#NOTES}

* The master branch includes the most stable version of GeN-Foam, but the develop branch include all the latest developments
* Please notice that models for water boiling are still preliminary, incomplete (missing models for boiling crisis) and in Beta testing.
* The adjoint diffusion solver has been implemented but not tested.
* The discrete-ordinate (SN) solver is currently only for steady-state. Transient analyses can in principle be run, but no acceleration techniques have been implemented, which requires hundreds of iterations for each time step.
* The removeBaffles flag in the controlDict may not always work in parallel. To be tested.
* The detailed temperature profile in the fuel cannot be recomposed with *reconstructPar*, as OpenFOAM is not currently able to reconstruct FieldFields. Be careful when trying to restart from a recomposed case, as you will have lost all info on fuel temperatures.
* In parallel cases, the keff in reactorState is updated only in the *processor* * subfolders.
