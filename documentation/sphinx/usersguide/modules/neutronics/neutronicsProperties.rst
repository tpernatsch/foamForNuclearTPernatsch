The *neutronicsProperties* dictionary
-------------------------------------

The *neutronicsProperties* dictionary is found under *constant/(nameOfNeutronicsRegion)/*
and is used to set the type of neutronics simulation by using the following keywords:

:eigenvalueNeutronics: should be set to ``true`` for eigenvalue calculations,
                       false for transients.
:externalSourceNeutronics: should be set to ``true`` for external neutron source
                           calculations. The ``eigenvalueNeutronics`` variable
                           should be put to ``false`` and the ``keff = 1``.
:fastNeutrons: should be set to ``true`` for fast spectrum reactors when using
               the ``pointKinetics`` solver.
:energyGroups: is an integer representing the number of energy groups used for
               3D spatial neutronics solver.
:precGroups: is an integer representing the number of delayed neutrons groups.
:adjustDiscFactors: when set to ``true``, the neutronics solver will try to
                    automatically adjust the discontinuity factors based on the
                    ``integralFlux`` provided in the ``nuclearData`` dictionary.
                    Please refer to Fiorina et al., Annals of Nuclear Energy
                    960 (2016)
:useGivenDiscFactors: when set to ``true``, the neutronics solver uses the
                      homogeneous discontinuity factors provided in the
                      ``nuclearData`` dictionary. Please refer to Fiorina et
                      al., Annals of Nuclear Energy 960 (2016)
:axialOrientation: fuel axial orientation vector used for thermomechanical
                   expansion.

One can find detailed, commented examples in most tutorials. See for instance
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/neutronicsProperties>`_ (single phase).

N.B: The parameter *model* used to define what type of simulation needs to be
performed has been replaced by the selection of model in the *system/regionsDict*
(see :ref:`Coupling solvers <userguide_coupling>`).
