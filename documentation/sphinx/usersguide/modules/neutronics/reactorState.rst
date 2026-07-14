.. _neutronics_the-reactorstate-dictionary:

The *reactorState* dictionary
-----------------------------

The *reactorState* dictionary is found under the *timeFolder/uniform/*
sub-folder. Its main keywords are:

:keff: is used in the spatial kinetics solvers as an initial guess for ``keff`` when
       doing an eigenvalue calculation. It is then updated automatically at each
       time step (i.e., at each power iteration) with the calculated value of
       ``keff``. When performing a transient calculation with the spatial kinetics
       solvers, ``keff`` is instead used to divide the neutron source term and is
       not updated during the simulation. Typically, to run spatial kinetics
       transient simulations, one first runs an eigenvalue calculation. The
       resulting ``keff`` will be the one that makes the reactor critical in a
       subsequent transient simulation. ``keff`` is disregarded by the point
       kinetics sub-solver.
:power: is used in the spatial kinetics solvers as target power when doing an
        eigenvalue calculation. It is also used by the point kinetics
        sub-solver, but only to correctly plot results. As power, the neutronics module
        uses what it finds under ``powerDensity`` of the *neutroRegion*, or under
        the ``powerDensity`` of the *fluidRegion* if it does not find a
        powerDensity in the neutroRegion. To correctly plot point kinetics
        results, ``power`` must be consistent with the mentioned power
        densities.
:precursorPowers: can be read by the point kinetics sub-solver in case the user
                  wishes to set initial concentrations of precursors. If not
                  found, precursor concentrations are initialized to be in
                  equilibrium with the starting conditions (i.e. a steady state
                  is assumed), :math:`C_i(0) = \frac{\beta_i}{\lambda_i \Lambda} \text{power}`.

The point kinetics sub-solver additionally writes and re-reads the reference
state used to evaluate the feedback reactivities, so that a restart resumes
from a consistent condition:

:TFuelRef, TCladRef, TCoolRef, rhoCoolRef, TStructRef, TStructMechRef, TDrivelineRef:
    the reference temperatures and coolant density against which the feedback
    reactivities are computed (:math:`\alpha_x (x - x^{ref})`). They are
    established at the first time step of a run (defaulting to the current
    field values if absent) and stored here. On a restart they are read back
    from *reactorState*, so the feedback baseline is preserved rather than
    reset to the restart-time values. Any of them may also be set by hand to
    impose a specific reference.
:precEquilibriumReactivity: (liquid fuel only) the reactivity offset that makes
    the circulating-fuel configuration critical at the start of the transient,
    :math:`\rho_{eq} = \beta_{tot} - \beta_{tot}^{eff}`. It is computed once,
    when the equilibrium precursor distribution is initialized
    (``initPrecursorsLiquidFuel true`` in *nuclearData*), and stored here. On a
    restart it is read back, so the (potentially hours-long) equilibrium solve
    does not have to be repeated: initialize once with the flag on, then restart
    subsequent transients with the flag off to reuse the stored value together
    with the persisted precursor fields. See the point-kinetics solver
    documentation for the full run-once / restart-many workflow.


All neutronics models can be used for liquid-fuel reactors. One can
activate this option using the ``liquidFuel`` keyword in */constant/(nameOfNeutronicsRegion)/neutronicsProperties*.
Of course, in such cases, one should pay attention to setting proper boundary
conditions for the precursors.

A commented *reactorState* can be found in `3D_SmallESFR
<https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/0/uniform/reactorState>`_.

.. note ::

    Please note that in parallel calculations, the updated *reactorState* can be 
    found in *timeStep/uniform/*.


.. warning ::

    The *powerDensity* file written to disk is the power density calculated
    by neutronics (``sigmaPow<g>`` multiplied by ``flux<g>``), DIVIDED by the ``fuelFraction``
    indicated in *nuclearData*. This means that it provides the power density IN the
    fuel, not spread over the cross-section homogenization region. For instance, if
    the model is an assembly with its own one-group cross-section set and the fuel fraction is specified
    0.3, ``powerDensity`` will be equal to:

    .. math::

        q''' = \frac{\sum_g \kappa_g \Sigma_{f,g} \phi_g}{\alpha_{fuel}} = \frac{\sum_g \text{sigmaPower}_g \times \text{flux}_g}{\text{fuelFraction}}
