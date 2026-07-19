The *neutronicsProperties* dictionary
-------------------------------------

The *neutronicsProperties* dictionary is found under *constant/(nameOfNeutronicsRegion)/*
and is used to set the type of neutronics simulation by using the following keywords:

:eigenvalueNeutronics: Flag to activate eigenvalue mode calculation.
                       If ``false`` performs transient calculation.
:externalSourceNeutronics: should be set to ``true`` for external neutron source
                           calculations. The ``eigenvalueNeutronics`` variable
                           should be put to ``false`` and the ``keff = 1``.
:fastNeutrons: Flag to activate reactivity feedback for fast spectrum reactors in 
               point-kinetics sub-solver. 
               If ``true``, use the ``-coeffDoppler * log(TFuel/TFuelRef)``, else 
               use ``coeffDoppler * (sqrt(TFuel) - sqrt(TFuelRef))``
:legendreMoments: Number of Legendre moments used for SN calculations.
:adjustDiscFactors: If ``true``, the neutronics sub-solver will try to
                    automatically adjust the discontinuity factors based on the
                    ``integralFlux`` provided in the ``nuclearData`` dictionary.
                    Please refer to :ref:`FIORINA2016212 <FIORINA2016212>`.
:useGivenDiscFactors: If ``true``, the neutronics sub-solver uses the
                      homogeneous discontinuity factors provided in the
                      ``nuclearData`` dictionary. Please refer to 
                      :ref:`FIORINA2016212 <FIORINA2016212>`.
:axialOrientation: Fuel axial orientation vector used for thermomechanical
                   expansion.
:groupsWoDF: List of energy groups that are not affected by discontinuity factors
             calculations.
:doNotParametrize: List of energy groups that are not parametrized.
:isLowMemory: Prepare low memory consumption for XS, useful for SN calculations.
:isReadXS: Is read XS flag.
:isWriteXS: Is write XS flag.


Here follows the subset of keyword for reactors with liquid fuel:

:liquidFuel: ``true`` for liquid fuel reactors calculations (e.g., Molten Salt Reactors),
             ``false`` for solid fuel.
:initPrecursorsLiquidFuel: Flag to initialize the concentration of precursors to their
                           equilibrium value at the beginning of the simulation (default 
                           ``false``).
:ScNo: Schmidt number for diffusion of precursors (default ``1``).


One can find detailed, commented examples in most tutorials. See for instance
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/neutronicsProperties>`_ (single phase).


.. note ::
 
    The parameter *model* used to define what type of simulation needs to be
    performed has been replaced by the selection of model in the *system/regionsDict*
    (see :ref:`Coupling solvers <couplingGF>`).
