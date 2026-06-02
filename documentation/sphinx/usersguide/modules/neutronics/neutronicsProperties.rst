The *neutronicsProperties* dictionary
-------------------------------------

The *neutronicsProperties* dictionary is found under *constant/(nameOfNeutronicsRegion)/*
and is used to set the type of neutronics simulation by using the following keywords:

:eigenvalueNeutronics: should be set to ``true`` for eigenvalue calculations,
                       false for transients.
:externalSourceNeutronics: should be set to ``true`` for external neutron source
                           calculations. The ``eigenvalueNeutronics`` variable
                           should be put to ``false`` and the ``keff = 1``.
:liquidFuel: should be set to ``true`` for liquid fuel reactors calculations (e.g., MSRs),
             false for solid fuel.
:fastNeutrons: should be set to ``true`` for fast spectrum reactors when using
               the ``pointKinetics`` solver.
:legendreMoments: Number of Legendre moments used for SN calculations.
:adjustDiscFactors: when set to ``true``, the neutronics solver will try to
                    automatically adjust the discontinuity factors based on the
                    ``integralFlux`` provided in the ``nuclearData`` dictionary.
                    Please refer to :ref:`FIORINA2016212 <FIORINA2016212>`.
:useGivenDiscFactors: when set to ``true``, the neutronics solver uses the
                      homogeneous discontinuity factors provided in the
                      ``nuclearData`` dictionary. Please refer to 
                      :ref:`FIORINA2016212 <FIORINA2016212>`.
:axialOrientation: fuel axial orientation vector used for thermomechanical
                   expansion.
:ScNo: Schmidt number.
:groupsWoDF: List of energy groups that are not affected by discontinuity factors
             calculations.
:doNotParametrize: List of energy groups that are not parametrized.
:isLowMemory: Prepare low memory consumption for XS, useful for SN calculations.
:isReadXS: Is read XS flag.
:isWriteXS: Is write XS flag.

One can find detailed, commented examples in most tutorials. See for instance
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/neutronicsProperties>`_ (single phase).


.. note ::
 
    The parameter *model* used to define what type of simulation needs to be
    performed has been replaced by the selection of model in the *system/regionsDict*
    (see :ref:`Coupling solvers <couplingGF>`).
