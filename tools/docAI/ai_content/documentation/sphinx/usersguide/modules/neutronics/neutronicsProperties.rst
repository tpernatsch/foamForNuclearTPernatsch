The *neutronicsProperties* dictionary
-------------------------------------

The *neutronicsProperties* dictionary is located under
``constant/(nameOfNeutronicsRegion)/``. It is used to set the type of
neutronics simulation using the following keywords:

:eigenvalueNeutronics: Flag to activate eigenvalue mode calculation.
                       If ``false``, perform a transient calculation.
:externalSourceNeutronics: Set to ``true`` for external neutron source
                           calculations. The ``eigenvalueNeutronics`` variable
                           should be set to ``false``, and ``keff`` should be set
                           to ``1``.
:fastNeutrons: Flag to activate reactivity feedback for fast spectrum
               reactors in the point-kinetics sub-solver.
               If ``true``, use ``-coeffDoppler * log(TFuel/TFuelRef)``.
               Otherwise, use ``coeffDoppler * (sqrt(TFuel) - sqrt(TFuelRef))``.
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
:groupsWoDF: List of energy groups that are not affected by discontinuity
             factor calculations.
:doNotParametrize: List of energy groups that are not parametrized.
:isLowMemory: Prepare low memory consumption for XS, useful for SN
               calculations.
:isReadXS: Flag to read XS.
:isWriteXS: Flag to write XS.


Here follows the subset of keywords for reactors with liquid fuel:

:liquidFuel: ``true`` for liquid fuel reactor calculations (e.g., Molten
             Salt Reactors). ``false`` for solid fuel.
:initPrecursorsLiquidFuel: Flag to initialize the concentration of precursors
                           to their equilibrium value at the beginning of the
                           simulation (default ``false``).
:ScNo: Schmidt number for diffusion of precursors (default ``1``).


One can find detailed, commented examples in most tutorials. See for instance
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/neutronicsProperties>`_ (single phase).


.. note ::

    The parameter *model* used to define what type of simulation needs to be
    performed has been replaced by the selection of model in the
    *system/regionsDict* (see :ref:`Coupling solvers <couplingGF>`).