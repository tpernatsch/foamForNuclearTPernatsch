==============
Setting a case
==============

Initial and boundary conditions
===============================

Initial and boundary conditions adopt the usual OpenFOAM logic for one- and
two-phase solvers. A couple of things to be kept in mind:

- The pressure field we solve for is *p_rgh* (pressure minus the gravitational head)
- When performing turbulent analyses, one needs to add the fields *k*, *epsilon*, nut and *alphat*

One thing that instead specific to GeN-Foam (except for the one-phase legacy
sub-solver) and that one needs to keep in mind is that U (or u.(name of fluid))
are the real velocities, not the Darcy velocities. In a porous structure, they
represent the actual velocity of the fluid, and not the velocity multiplied by
the fluid fraction. For instance, U will increase when transiting from a high-
to a low-porosity region.

OpenFOAM provides most of the boundary conditions one may need for
thermal-hydraulics models. In addition, a few boundary conditions have been
included in GeN-Foam and can be found in
`GeN-Foam/classes/thermalHydraulics/src/boundaryConditions <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/src/fvPatchFields/thermalHydraulics>`_.
Information on the use of each boundary condition can be found in the header
files (.H).


.. _userguide_thermalhydraulics_settingcase_setpower:

Setting the initial power
-------------------------

There are several ways to set the power in GeN-Foam. For the power generated in
subscale structures:

- The thermal-hydraulics solver will normally use the *powerDensity.* fields (for instance, *powerDensity.nuclearFuelPin* for pin-based reactors) that it finds in the 0 (or *startTime*) folder.
- As an alternative, one can provide cellZone-by-cellZone values via the keyword *powerDensity* in the various power models in the *phase* properties (see for instance `1D_boiling <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/1D_boiling/constant/fluidRegion/phaseProperties>`_). However, if GeN-Foam finds the corresponding field in the 0 (or *startTime*) folder, this will take priority.

For the power generated in the fluid itself:

- The thermal-hydraulics solver will normally use the *powerDensity* field that it finds in the 0 (or *startTime*) folder.
- One can override this behavior by using the *initialPowerDensity* keyword in the *phaseProperties* (see `1D_MSR_pointKinetics <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/1D_MSR_pointKinetics/rootCase/constant/fluidRegion/phaseProperties>`_) for an example. Also in this case the field in the 0 (or *startTime*) folder will take priority.

If neutronics is activated, the power density can be mapped from the neutronics
sub-solver. For eigenvalue calculations, the power is set in the *pTarget*
keyword in the *reactorState* dictionary. For transients, the power is a result
of calculations. There is one important exception to this behavior: the point
kinetics solver will only rescale the power densities (see below about why
plural) that it finds in *neutroRegion*. Alternatively, the powerDensity field
can be constructed by the thermal-hydraulics solver and mapped once from
*fluidRegion* to *neutroRegion*. It is possible to achieve this either using the
*mapFields* utility of OpenFOAM or by running
``GeN-Foam -initializeMappedFields``, which will construct the solvers and
mappings based on the *multiRegionCouplingDict* and map the fields without
running the simulation. For point kinetics, the *pTarget* keyword in
*reactorState* is not used by the solver itself. However, to correctly plot
point kinetics results, pTarget must be consistent with the mentioned power
densities.

.. note::

    In two-phase simulations with liquid fuel, the powerDensity in neutronics
    goes to anything that is liquid in thermal-hydraulics. You are supposed to
    have one liquid and one gas. Otherwise, power will be counted twice.


Power densities and secondary power densities, and liquid fuel
--------------------------------------------------------------

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
===========================

Details for discretization and solution of equations are handled in a standard
OpenFOAM way, i.e., through the *fvSolution* and *fvSchemes* dictionaries in
*system/fluidRegion*.
