The *neutronicsProperties* dictionary
-------------------------------------

The *neutronicsProperties* dictionary is found under *constant/(nameOfNeutronicsRegion)/*
and is used to set the type of neutronics simulation by using the following keywords:

:eigenvalueNeutronics: should be set to ``true`` for eigenvalue calculations,
                       false for transients.
:externalSourceNeutronics: should be set to ``true`` for external neutron source
                           calculations. The ``eigenvalueNeutronics`` variable
                           should be put to ``false`` and the ``keff = 1``.

One can find detailed, commented examples in most tutorials. See for instance
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/neutronicsProperties>`_ (single phase).

N.B: The parameter *model* used to define what type of simulation needs to be
performed as been replaced by the selection of model in the *system/regionsDict*
(see :ref:`Coupling solvers <userguide_coupling>`).


Models
~~~~~~

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


