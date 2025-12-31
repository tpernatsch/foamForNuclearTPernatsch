.. _neutronics_tipsAndTricks:



.. note ::
    ``defaultPrec`` field has 1/m3 units except for the adjoint solver that needs 1/m2/s.

.. note ::
    In point kinetics, pTarget (in reactorState) MUST be the same as the one used for reaching the steady-state. As power, GeN-Foam uses what it finds under powerDensity, or under the powerDensity of the fluidRegion if it does not find a powerDensity in the neutroRegion. pTarget does not enter the calculation, it is used simply to plot the results