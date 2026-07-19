.. _userguide_neutronics_subsolvers:

Sub-solvers
-----------

Neutronics calculations are performed by classes derived from
:ref:`neutronics.H <neutronicsModel>` that contain specific sub-solvers.
For the user, the derived classes translate into runtime-selectable models.
The specific sub-solver used in a simulation is normally selected at the level
of the application that uses the modules. For example, in GeN-Foam, the module
is chosen in the ``regionsDict`` (see :ref:`Achieving coupled solutions <couplingGF>`).

Neutronics calculations are performed by classes derived from *neutronics* that
contain specific sub-solvers:

.. toctree::
   :maxdepth: 1
   :glob:

   ../../../cppapi/generated/modules/neutronics/**


.. raw:: html

   <br><br>

``adjointDiffusion`` has been developed only as an eigenvalue solver. The others
can be used for transient calculations. However, the SN transient solver has not
been tested. In addition, it is currently not accelerated, and is therefore
extremely slow. It can require hundreds of iterations per time step.


Subcritical point-kinetics
~~~~~~~~~~~~~~~~~~~~~~~~~~

To use the subcritical point-kinetics, the user must add
*constant/neutroRegion/externalSource*. The file contains a flag to activate the
external neutron source (``isExternalSource``).

Several parameters related to a spallation source are included, such as the
energy per source particle (in J/source particle) and the neutron yield of the
reaction (in neutrons/source particle).

An external source modulation timetable is provided to manually modulate the
source strength.

In the case of an FMI coupling, it is possible to use the
``externalSourceModulationNameFromFMU`` entry to change the external source
modulation through an FMI. To use it, the mode must be ``transient``.


Power normalization
~~~~~~~~~~~~~~~~~~~

The spatial neutronics solvers always create the *powerDensity* and
*secondaryPowerDensity* fields. By default, *secondaryPowerDensity* is set to
zero, and the ``fuelFraction`` keyword in *nuclearData* is used to translate the
volume-average power density (normally calculated by multiplying cross-sections
and fluxes) into the fuel-averaged power density required by the
thermal-hydraulic sub-solver.

However, a *secondaryPowerDensity* may sometimes be needed. It can be used
to provide power to the coolant in a solid-fuel reactor and, more importantly,
to provide power to the graphite in a liquid-fuel reactor.

To calculate a *secondaryPowerDensity*, the neutronics solver needs to know
how much of the total power goes into *secondaryPowerDensity*, and the volume
fraction of the secondary power-producing structure or liquid. This can be done
by using the ``fractionToSecondaryPower`` and ``secondaryPowerVolumeFraction``
keywords in each cellZone in the *nuclearData* sub-dictionary (the same place as
``fuelFraction``). If these keywords are present, the neutronics module will
calculate power densities as follows:

- :math:`\text{secondaryPowerDensity} = \frac{\text{powerDensity}}{\max(\text{secondaryPowerVolumeFraction}, \text{SMALL})} \times \text{fractionToSecondaryPower}`
- :math:`\text{powerDensity} = \frac{\text{powerDensity}}{\max(\text{fuelFraction}, SMALL)} \times (1.0 - \text{fractionToSecondaryPower})`


When the point kinetics sub-solver is selected, the *powerDensity* and
*secondaryPowerDensity* fields found in the initial time folder are rescaled.
The only exception is when ``liquidFuel`` is true and the ``initialPowerDensity``
keyword is used. In this case, ``initialPowerDensity`` takes priority. This is
the value that the neutronics sub-solver will rescale and print.