.. _neutronics_the-reactorstate-dictionary:

The *reactorState* dictionary
-----------------------------

The *reactorState* dictionary is found under the *timeFolder/uniform/*
sub-folder. It includes essentially 3 keywords:

:keff: is used in the spatial kinetics solvers as an initial guess for *keff* when
       doing an eigenvalue calculation. It is then updated automatically at each
       time step (i.e., at each power iteration) with the calculated value of
       *keff*. When performing a transient calculation with the spatial kinetics
       solvers, *keff* is instead used to divide the neutron source term and is
       not updated during the simulation. Typically, to run spatial kinetics
       transient simulations, one first runs an eigenvalue calculation. The
       resulting *keff* will be the one that makes the reactor critical in a
       subsequent transient simulation. *keff* is disregarded by the point
       kinetics sub-solver.
:power: is used in the spatial kinetics solvers as target power when doing an
        eigenvalue calculation. It is also used by the point kinetics
        sub-solver, but only to correctly plot results. As power, GeN-Foam
        uses what it finds under powerDensity of the neutroRegion, or under
        the powerDensity of the fluidRegion if it does not find a
        powerDensity in the neutroRegion. To correctly plot point kinetics
        results, pTarget must be consistent with the mentioned power
        densities.
:precursorPowers: can be read by the point kinetics sub-solver in case the user
                  wishes to set initial concentrations of precursors. If not
                  found, precursor concentrations are initialized to be in
                  equilibrium with the starting conditions (i.e. a steady state
                  is assumed).


All neutronics models can be used for liquid-fuel reactors. One can
activate this option using the *liquidFuel* keyword in */system/controlDict*.
Of course, in such cases, one should pay attention to setting proper boundary
conditions for the precursors.

A commented *reactorState* can be found in `3D_SmallESFR
<https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/0/uniform/reactorState>`_.

.. note ::
    Please note that in parallel calculations, the updated *reactorState* can be 
    found in *timeStep/uniform/*.


.. warning ::
    The powerDensity file written to disk is the power density calculated
    by neutronics (sigmaPowers multiplied by fluxes), DIVIDED by the fuelFractions
    indicated in nuclearData. This means that it provides the power density IN the
    fuel, not spread over the cross-section homogenization region. For instance, if
    you have an assembly with its own one-group cross-section set and you specify
    that the fuel fraction is 0.3, powerDensity will be equal to:

    .. math::

        q''' = \frac{\kappa \Sigma_f \phi}{\alpha_{fuel}} = \frac{\text{sigmaPower} \times \text{flux}}{\text{fuelFraction}}
