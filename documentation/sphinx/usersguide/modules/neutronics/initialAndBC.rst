.. _userguide_thermalhydraulics_initialAndBC:


Initial and boundary conditions
-------------------------------

Initial and boundary conditions adopt the usual OpenFOAM logic for one- and two-phase solvers. 
OpenFOAM provides most of the boundary conditions one may need for
thermal-hydraulics models. In addition, a few boundary conditions have been
included in the thermal-hydraulics module of foamForNuclear:

toctreeHere


For the initial conditions of fluxes, 
the user can either specify different IC and BC for each one of the energy
groups (with fluxes that must be named *fluxStar0*, *fluxStar1*, etc…) or
provide the same IC and BC to all fluxes by using the *defaultFlux* field. In
the case of SP3 calculations, the IC and BC for the second moment can be imposed
either for each energy (using fields named *fluxStar20*, *fluxStar21*, etc…),
or to all energies by using the *defaultFlux2* field. When both *defaultFlux*
and *fluxStar...* are present, the solver gives priority to *fluxStar...*. In
the case of SN calculations, it is suggested not to modify the boundary
conditions and to use the *defaultFlux* file (an example is provided in the
Godiva_SN tutorial). When employing the adjoint solver, the user will have to add the
fields *adjointDefaultPrec* and *adjointDefaultFlux* in the  initial time.


IC and BC for precursors do not have to be specified for standard reactors. On
the other hand, they should be specified in the case of liquid fuel reactors
(e.g., Molten Salt Reactors). This is possible by creating a *defaultPrec*
field, in case the same conditions apply to all precursor groups, or by creating
the fields named *prec0*, *prec1*, etc., in case different conditions must be
provided for different precursor groups.

.. note ::

    Boundary conditions must be applied to ``fluxStar...`` and not to ``flux...``
    since GeN-Foam solves for these variables. ``fluxStar...`` represent continuous
    fluxes, while ``flux...`` represent the real fluxes. They differ only in case
    discontinuity factors are employed (see :ref:`FIORINA2016212 <FIORINA2016212>`).


Setting the weighting in point-kinetics calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A correct evaluation of the reactivity worth of delayed neutron precursors in
MSRs, as well as of the impact of temperatures on reactivities, normally
requires the knowledge of the adjoint flux. In GeN-Foam, the field
*oneGroupFlux* is used by the point kinetic solver for weighting temperatures,
densities and precursors. When fluxes are not calculated via a spatial
neutronics calculation, one has to manually provide the *oneGroupFlux* in
*0/neutroRegion*. As an alternative, one can use the *initialOneGroupFluxByZone*
keyword in *nuclearData* (see `1D_MSR_pointKinetics
<https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/1D_MSR_pointKinetics/rootCase/constant/neutroRegion/nuclearData>`_). Please notice that:

- If calculated fluxes are available in *neutroRegion*, these will be user to recalculate and overwrite *oneGroupFlux*.
- If no fluxes are available, the neutronics sub-solver will use the provided *oneGroupFlux*
- If the *initialOneGroupFluxByZone* keyword is used in *nuclearData*, this will be used to overwrite *oneGroupFlux*


Power densities and secondary power densities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
