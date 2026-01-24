.. _userguide_thermalhydraulics_settingcase:

Setting a case
--------------

Initial and boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Initial and boundary conditions adopt the usual OpenFOAM logic for one- and two-phase solvers. 
OpenFOAM provides most of the boundary conditions one may need for
thermal-hydraulics models. In addition, a few boundary conditions have been
included in GeN-Foam and can be found in
`GeN-Foam/classes/thermalHydraulics/src/boundaryConditions <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/src/fvPatchFields/thermalHydraulics>`_.
Information on the use of each boundary condition can be found in the header
files (.H).




.. note::

    In two-phase simulations with liquid fuel, the powerDensity in neutronics
    goes to anything that is liquid in thermal-hydraulics. You are supposed to
    have one liquid and one gas. Otherwise, power will be counted twice.


Power densities and secondary power densities, and liquid fuel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The spatial neutronics solvers always create a *powerDensity* and
*secondaryPowerDensity* fields. By default, *secondaryPowerDensity* is set to
zero and the *fuelFraction* keyword in *nuclearData* is used to translate the
volume-average power density that is normally calculated by multiplying
cross-sections and fluxes into the fuel-averaged power density that is needed by
the thermal-hydraulic sub-solver

However, a *secondaryPowerDensity* might sometimes be needed. It might be used
to provide some power to the coolant in a solid-fuel reactor and, more
importantly, to provide some power to the graphite in a liquid-fuel reactor. To
calculate a *secondaryPowerDensity*, GeN-Foam needs to know how much of the
total power goes into the *secondaryPowerDensity*, and what is the volume
fraction of the secondary power-producing structure or liquid. This can be done
by using the *fractionToSecondaryPower* and *secondaryPowerVolumeFraction*
keywords in each cellZone in *nuclearData* (the same place as *fuelFraction*).
If these keywords are present, GeN-Foam will calculate power densities as
follows:

- :math:`\text{secondaryPowerDensity} = \frac{\text{powerDensity}}{\max(\text{secondaryPowerVolumeFraction}, \text{SMALL}) \times \text{fractionToSecondaryPower}}`
- :math:`\text{powerDensity} = \frac{\text{powerDensity}}{\max(\text{fuelFraction}, SMALL) \times (1.0 - \text{fractionToSecondaryPower})}`

When the *liquidFuel* flag is set to false, the thermal-hydraulic sub-solver will:

- take the *powerDensity* field from neutronics and project it to its own *powerDensityNeutronics* field;
- take the *secondaryPowerDensity* field from neutronics and project it to its own *powerDensityNeutronicsToLiquid* field.

When the *liquidFuel* flag is set to true, the thermal-hydraulic sub-solver will:

- take the *powerDensity* field from neutronics and project it to its own *powerDensityNeutronicsToLiquid* field;
- take the *secondaryPowerDensity* field from neutronics and project it to its own *powerDensityNeutronics* field.

When point kinetics is used, the solver will simply rescale the *powerDensity*
and *secondaryPowerDensity* it finds, and the thermal-hydraulic solver will take
them depending on the *liquidFuel* flag as described above. The only exception
is when *liquidFuel* is true and the *initialPowerDensity* keyword is used. In
this case, *initialPowerDensity* will take priority and this is the value that
GeN-Foam will rescale and print to the powerDensityToLiquid

.. note::

    The power densities in the thermal-hydraulic sub-solver are ALWAYS the
    physical ones: for instance, when the *nuclearFuelPin* model is used for
    pin-based reactors, *powerDensity* refers to the power density inside the
    fuel matrix. For liquid fuel, the *powerDensity* is the power density in the
    liquid. They are not the power densities smeared over the whole volume.

.. note::

    If you calculate the powerDensity using Serpent, you have to divide it by
    the fuel fraction before feeding it to GeN-Foam (see
    :ref:`Neutronics <neutronics_the-reactorstate-dictionary>`)


Discretization and solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Details for discretization and solution of equations are handled in a standard
OpenFOAM way, i.e., through the *fvSolution* and *fvSchemes* dictionaries in
*system/fluidRegion*.
