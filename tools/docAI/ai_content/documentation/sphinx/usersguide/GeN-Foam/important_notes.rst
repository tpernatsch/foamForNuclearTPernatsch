.. _importantnotes:


Important notes
===============

- The ``master`` branch includes the most stable version of GeN-Foam, while the ``develop`` branch includes the latest developments.
- Please note that the models for water boiling are still preliminary, incomplete, and in beta testing. In particular, the models for CHF and post-CHF have not been verified yet. Use with care. For the moment, the critical heat flux can only be imposed as a constant value. No models or lookup tables have been implemented yet for its prediction.
- The adjoint diffusion solver has been implemented, but it has not been tested.
- The discrete-ordinate (SN) solver is currently limited to steady-state calculations. Transient analyses can, in principle, be run, but no acceleration techniques have been implemented. This requires hundreds of iterations for each time step.
- The ``removeBaffles`` flag in ``controlDict`` may not always work in parallel. This needs to be tested.
- The detailed temperature profile in the fuel cannot be recomposed with ``reconstructPar``, because OpenFOAM is currently not able to reconstruct ``FieldFields``. Be careful when trying to restart from a recomposed case, as you will have lost all information on fuel temperatures.