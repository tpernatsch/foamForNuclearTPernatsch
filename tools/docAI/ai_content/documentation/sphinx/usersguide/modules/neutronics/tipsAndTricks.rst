.. _neutronics_tipsAndTricks:

Tips and tricks
---------------

.. note ::

    The ``defaultPrec`` field has units of 1/m3, except for the adjoint solver,
    which requires units of 1/m2/s.

.. note ::

    In point kinetics, ``power`` (in *reactorState*) MUST be the same as the
    one used to reach the steady state. For power, the neutronics module uses
    what it finds under *powerDensity*, or under the *powerDensity* of the
    *fluidRegion* if it does not find a *powerDensity* in the *neutroRegion*.
    ``power`` does not enter the calculation. It is used only to plot the
    results.