.. _usersguide_modules_neutronics:


Neutronics
==========

The neutronics module  (see :ref:`neutronics.H <neutronics>`)  can be used to
solve for steady-state and transient neutronics according to different
approximations. The behavior of the neutronics module is dictated by the
*neutronicsProperties*,  *nuclearData*, *reactorState*, *quadratureSet*, and
*quadratureSet* dictionaries. 

Note that the spatial neutronics solvers always create a *powerDensity* and
*secondaryPowerDensity* fields. By default, *secondaryPowerDensity* is set to
zero and the *fuelFraction* keyword in *nuclearData* is used to translate the
volume-average power density that is normally calculated by multiplying
cross-sections and fluxes into the fuel-averaged power density that is needed by
the thermal-hydraulic sub-solver. However, a *secondaryPowerDensity* might sometimes be needed. It might be used
to provide some power to the coolant in a solid-fuel reactor and, more
importantly, to provide some power to the graphite in a liquid-fuel reactor. To
calculate a *secondaryPowerDensity*, GeN-Foam needs to know how much of the
total power goes into the *secondaryPowerDensity*, and what is the volume
fraction of the secondary power-producing structure or liquid. This can be done
by using the *fractionToSecondaryPower* and *secondaryPowerVolumeFraction*
keywords in each cellZone in the *nuclearData* sub-dictionaty (the same place as *fuelFraction*).
If these keywords are present, GeN-Foam will calculate power densities as
follows:

- :math:`\text{secondaryPowerDensity} = \frac{\text{powerDensity}}{\max(\text{secondaryPowerVolumeFraction}, \text{SMALL})} \times \text{fractionToSecondaryPower}`
- :math:`\text{powerDensity} = \frac{\text{powerDensity}}{\max(\text{fuelFraction}, SMALL)} \times (1.0 - \text{fractionToSecondaryPower})`


When point kinetics is used, the solver will simply rescale the *powerDensity*
and *secondaryPowerDensity* it finds in the initial time folder. The only exception
is when *liquidFuel* is true and the *initialPowerDensity* keyword is used. In
this case, *initialPowerDensity* will take priority and this is the value that
GeN-Foam will rescale and print.


.. toctree::
    :maxdepth: 1

    neutronicsSolvers
    initialAndBC
    discretizationSolution
    neutronicsProperties
    nuclearData
    reactorState
    quadratureSet
    CRMove
    tipsAndTricks
