.. _neutronics_the-reactorstate-dictionary:

The *reactorState* dictionary
-----------------------------

The *reactorState* dictionary is found under the *timeFolder/uniform/*
sub-folder. It includes essentially 3 keywords:

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
