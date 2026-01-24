
Sub-solvers
-----------

Neutronics calculations are performed by classes derived from *neutronics* that
contain specific sub-solvers:

:pointKinetics: Point-kinetics (:ref:`pointKineticNeutronics.H <pointKineticNeutronics>`)
:diffusionNeutronics: Diffusion (:ref:`diffusionNeutronics.H <diffusionNeutronics>`)
:adjointDiffusion: Adjoint diffusion (:ref:`adjointDiffusionNeutronics.H <adjointDiffusionNeutronics>`)
:SP3Neutronics: Diffusion in :math:`SP_3` (:ref:`SP3Neutronics.H <SP3Neutronics>`)
:SNNeutronics: Discrete ordinates (:math:`S_N`) (:ref:`SNNeutronics.H <SNNeutronics>`)

For the user, the derived classes translate into runtime selectable models. The
specific sub-solver to be used in a simulation can be selected at runtime in the
*constant/neutroRegion/neutronicsProperties* dictionary.

The choice of the model is achieved by selecting the wanted solver in
*regionSolvers* depending on whether the neutronics
solvers need to be part of a tightly coupled loop or not (see
:ref:`Coupling solvers <userguide_coupling_the-controlDict-dictionary>`).

*adjointDiffusion* has been developed only as an eigenvalue solver. The others
can be used for transient calculations. However, the SN transient solver has not
been tested. In addition, it is currently not accelerated, thus extremely slow
(it can require hundreds of iterations per time step).


Subcritical point-kinetics
~~~~~~~~~~~~~~~~~~~~~~~~~~

To use the subcritical point-kinetics, the user has to add the
*constant/neutroRegion/externalSource*. The file contains a flag to activate the
external neutron source (``isExternalSource``).

Several parameters related to a spallation source are included such as the
energy per source particle in J/source particle and the neutron yield of the
reaction in neutrons/source particle*

An external source modulation timetable is provided to manually modulate the
source strength.

In the case of an FMI coupling, it is possible to use the
``externalSourceModulationNameFromFMU`` entry to change the external source
modulation through an FMI. To use it, the mode must be ``transient``.