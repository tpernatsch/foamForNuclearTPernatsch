.. _cppapi:

===============
C++ Source Code
===============

----------
Neutronics
----------

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/modules/neutronics/neutronics
   generated/modules/neutronics/SP3/SP3Neutronics
   generated/modules/neutronics/diffusion/diffusionNeutronics
   generated/modules/neutronics/adjointDiffusion/adjointDiffusionNeutronics
   generated/modules/neutronics/pointKinetics/pointKineticNeutronics
   generated/modules/neutronics/SN/SNNeutronics

.. list-table::
    :widths: 50 50

    * - :doc:`generated/modules/neutronics/neutronics`
      - Base class for neutronic models
    * - :doc:`generated/modules/neutronics/SP3/SP3Neutronics`
      - Subclass for SP3 neutronics
    * - :doc:`generated/modules/neutronics/diffusion/diffusionNeutronics`
      - Subclass for diffusion neutronics
    * - :doc:`generated/modules/neutronics/adjointDiffusion/adjointDiffusionNeutronics`
      - Derived class for adjointDiffusion neutronics
    * - :doc:`generated/modules/neutronics/pointKinetics/pointKineticNeutronics`
      - Run-time selectable point-kinetic based neutronics solver :ref:`RADMAN2022108891 <RADMAN2022108891>`
    * - :doc:`generated/modules/neutronics/SN/SNNeutronics`
      - Subclass for SN neutronics (Discrete Ordinate)


------------------
Thermal-hydraulics
------------------



1-phase Solvers
---------------

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/modules/thermalHydraulics/onePhaseRhoPimpleFoam/rhoPimpleFoam
   generated/modules/thermalHydraulics/onePhase/onePhase
   generated/modules/thermalHydraulics/onePhaseLegacy/onePhaseLegacy

.. list-table::
    :widths: 50 50

    * - :doc:`generated/modules/thermalHydraulics/onePhaseRhoPimpleFoam/rhoPimpleFoam`
      - This class encapsulates a the OpenFOAM solver rhoPimpleFoam
    * - :doc:`generated/modules/thermalHydraulics/onePhase/onePhase`
      - This class encapsulates a single fluid + 1 stationary structure
    * - :doc:`generated/modules/thermalHydraulics/onePhaseLegacy/onePhaseLegacy`
      - This class is based on onePhase


2-phase Solvers
---------------

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/modules/thermalHydraulics/twoPhase/twoPhase

.. list-table::
    :widths: 50 50

    * - :doc:`generated/modules/thermalHydraulics/twoPhase/twoPhase`
      - This class encapsulates a two fluid + 1 stationary structure


Phase Models
------------



Power Models
^^^^^^^^^^^^

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/porousMediaModels/phaseModels/structureModels/powerModels/powerModel
   generated/porousMediaModels/phaseModels/structureModels/powerModels/fixedPower/fixedPower
   generated/porousMediaModels/phaseModels/structureModels/powerModels/fixedTemperature/fixedTemperature
   generated/porousMediaModels/phaseModels/structureModels/powerModels/nuclearFuelPin/nuclearFuelPin
   generated/porousMediaModels/phaseModels/structureModels/powerModels/interpolatedNuclearFuelPin/interpolatedNuclearFuelPin
   generated/porousMediaModels/phaseModels/structureModels/powerModels/heatedPin/heatedPin
   generated/porousMediaModels/phaseModels/structureModels/powerModels/nuclearFuelFMU/nuclearFuelFMU
   generated/porousMediaModels/phaseModels/structureModels/powerModels/nuclearFuelFMU/listOperation/listConversion
   generated/porousMediaModels/phaseModels/structureModels/powerModels/nuclearFuelFMU/latticeMap/latticeMap
   generated/porousMediaModels/phaseModels/structureModels/powerModels/nuclearFuelPinTest/nuclearFuelPinTest
   generated/porousMediaModels/phaseModels/structureModels/powerModels/fixedTemperatureFMU/fixedTemperatureFMU
   generated/porousMediaModels/phaseModels/structureModels/powerModels/lumpedNuclearStructure/lumpedNuclearStructure
   generated/porousMediaModels/phaseModels/structureModels/powerModels/nuclearSteadyStatePebble/nuclearSteadyStatePebble

.. list-table::
    :widths: 50 50

    * - :doc:`generated/porousMediaModels/phaseModels/structureModels/powerModels/powerModel`
      - Run-time selectable class to handle the thermal description (i
    * - :doc:`generated/porousMediaModels/phaseModels/structureModels/powerModels/fixedPower/fixedPower`
      - Describes a generic structure with a constant or time depdendent internal
    * - :doc:`generated/porousMediaModels/phaseModels/structureModels/powerModels/fixedTemperature/fixedTemperature`
      - Describes a structure with a constant surface temperature
    * - :doc:`generated/porousMediaModels/phaseModels/structureModels/powerModels/nuclearFuelPin/nuclearFuelPin`
      - Model for representing a nuclear fuel pin with constant material
    * - :doc:`generated/porousMediaModels/phaseModels/structureModels/powerModels/interpolatedNuclearFuelPin/interpolatedNuclearFuelPin`
      - Model for representing a nuclear fuel pin with constant material
    * - :doc:`generated/porousMediaModels/phaseModels/structureModels/powerModels/heatedPin/heatedPin`
      - Model for representing a heated pin with constant material properties that
    * - :doc:`generated/porousMediaModels/phaseModels/structureModels/powerModels/nuclearFuelFMU/nuclearFuelFMU`
      - Model for representing a nuclear fuel from an external FMU
    * - :doc:`generated/porousMediaModels/phaseModels/structureModels/powerModels/nuclearFuelFMU/listOperation/listConversion`
      - Base file for list conversion
    * - :doc:`generated/porousMediaModels/phaseModels/structureModels/powerModels/nuclearFuelFMU/latticeMap/latticeMap`
      - Base class for lattices
    * - :doc:`generated/porousMediaModels/phaseModels/structureModels/powerModels/nuclearFuelPinTest/nuclearFuelPinTest`
      - Same as nuclearFuelPin
    * - :doc:`generated/porousMediaModels/phaseModels/structureModels/powerModels/fixedTemperatureFMU/fixedTemperatureFMU`
      - Describes a structure with a temperature taken from an FMU
    * - :doc:`generated/porousMediaModels/phaseModels/structureModels/powerModels/lumpedNuclearStructure/lumpedNuclearStructure`
      - Model for representing a lumped-parameter structure characterized
    * - :doc:`generated/porousMediaModels/phaseModels/structureModels/powerModels/nuclearSteadyStatePebble/nuclearSteadyStatePebble`
      - Describes triso and pebble of a pebble-bed reactor


Heat Exchanger
^^^^^^^^^^^^^^

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/porousMediaModels/phaseModels/structureModels/heatExchanger/heatExchanger

.. list-table::
    :widths: 50 50

    * - :doc:`generated/porousMediaModels/phaseModels/structureModels/heatExchanger/heatExchanger`
      - Heat exchanger


Pump
^^^^

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/porousMediaModels/phaseModels/structureModels/pump/pump

.. list-table::
    :widths: 50 50

    * - :doc:`generated/porousMediaModels/phaseModels/structureModels/pump/pump`
      - Allow defining a momentum source, including options for a time


Power Off Criterion
^^^^^^^^^^^^^^^^^^^

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/porousMediaModels/phaseModels/structureModels/powerOffCriterionModels/powerOffCriterionModel
   generated/porousMediaModels/phaseModels/structureModels/powerOffCriterionModels/fieldValue/powerOffCriterionFieldValue
   generated/porousMediaModels/phaseModels/structureModels/powerOffCriterionModels/timer/powerOffCriterionTimer

.. list-table::
    :widths: 50 50

    * - :doc:`generated/porousMediaModels/phaseModels/structureModels/powerOffCriterionModels/powerOffCriterionModel`
      - Author:
    * - :doc:`generated/porousMediaModels/phaseModels/structureModels/powerOffCriterionModels/fieldValue/powerOffCriterionFieldValue`
      - Turns power off if the field value maximum exceeds a provided value
    * - :doc:`generated/porousMediaModels/phaseModels/structureModels/powerOffCriterionModels/timer/powerOffCriterionTimer`
      - Turns power off if the field value maximum exceeds a provided value


Physics Models
--------------



Drag Models
^^^^^^^^^^^

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragFactor
   generated/porousMediaModels/physicsModels/dragModels/FFDragFactor
   generated/porousMediaModels/physicsModels/dragModels/FSDragFactor
   generated/porousMediaModels/physicsModels/dragModels/FFDragCoefficientModels/FFDragCoefficientModel
   generated/porousMediaModels/physicsModels/dragModels/FFDragCoefficientModels/Autruffe/AutruffeFFDragCoefficient
   generated/porousMediaModels/physicsModels/dragModels/FFDragCoefficientModels/NoKazimi/NoKazimiFFDragCoefficient
   generated/porousMediaModels/physicsModels/dragModels/FFDragCoefficientModels/BestionTRACE/BestionTRACEFFDragCoefficient
   generated/porousMediaModels/physicsModels/dragModels/FFDragCoefficientModels/Wallis/WallisFFDragCoefficient
   generated/porousMediaModels/physicsModels/dragModels/FFDragCoefficientModels/SchillerNaumann/SchillerNaumannFFDragCoefficient
   generated/porousMediaModels/physicsModels/dragModels/FFDragCoefficientModels/Bestion/BestionFFDragCoefficient
   generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragMultiplierModels/twoPhaseDragMultiplierModel
   generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragMultiplierModels/Kaiser88/Kaiser88TwoPhaseDragMultiplier
   generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragMultiplierModels/KottowskiSavatteri/KottowskiSavatteriTwoPhaseDragMultiplier
   generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragMultiplierModels/LockhartMartinelli/LockhartMartinelli
   generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragMultiplierModels/ChenKalish/ChenKalishTwoPhaseDragMultiplier
   generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragMultiplierModels/LottesFlinn/LottesFlinnTwoPhaseDragMultiplier
   generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragMultiplierModels/LottesFlinnNguyen/LottesFlinnNguyenTwoPhaseDragMultiplier
   generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragMultiplierModels/constant/constantTwoPhaseDragMultiplier
   generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragMultiplierModels/Kaiser74/Kaiser74TwoPhaseDragMultiplier
   generated/porousMediaModels/physicsModels/dragModels/FSDragCoefficientModels/FSDragCoefficientModel
   generated/porousMediaModels/physicsModels/dragModels/FSDragCoefficientModels/Rehme/RehmeFSDragCoefficient
   generated/porousMediaModels/physicsModels/dragModels/FSDragCoefficientModels/Churchill/ChurchillFSDragCoefficient
   generated/porousMediaModels/physicsModels/dragModels/FSDragCoefficientModels/Colebrook/ColebrookFSDragCoefficient
   generated/porousMediaModels/physicsModels/dragModels/FSDragCoefficientModels/NoKazimi/NoKazimiFSDragCoefficient
   generated/porousMediaModels/physicsModels/dragModels/FSDragCoefficientModels/modifiedEngel/modifiedEngelFSDragCoefficient
   generated/porousMediaModels/physicsModels/dragModels/FSDragCoefficientModels/BaxiDalleDonne/BaxiDalleDonneFSDragCoefficient
   generated/porousMediaModels/physicsModels/dragModels/FSDragCoefficientModels/Engel/EngelFSDragCoefficient
   generated/porousMediaModels/physicsModels/dragModels/FSDragCoefficientModels/ReynoldsPower/ReynoldsPowerFSDragCoefficient

.. list-table::
    :widths: 50 50

    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragFactor`
      - Author:
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/FFDragFactor`
      - Base class to represent a FFDragFactor with basing name access
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/FSDragFactor`
      - Base class to represent a FSDragFactor with basing name access
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/FFDragCoefficientModels/FFDragCoefficientModel`
      - Class to manage the fluid-structure (i
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/FFDragCoefficientModels/Autruffe/AutruffeFFDragCoefficient`
      - Autruffe correlation for the friction factor between two fluids,
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/FFDragCoefficientModels/NoKazimi/NoKazimiFFDragCoefficient`
      - No-Kazimi model for interfacial friction
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/FFDragCoefficientModels/BestionTRACE/BestionTRACEFFDragCoefficient`
      - Bestion model for interfacial friction
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/FFDragCoefficientModels/Wallis/WallisFFDragCoefficient`
      - Wallis correlation for interfacial liquid-vapour sodium friction
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/FFDragCoefficientModels/SchillerNaumann/SchillerNaumannFFDragCoefficient`
      - Schiller-Neumann model for interfacial friction
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/FFDragCoefficientModels/Bestion/BestionFFDragCoefficient`
      - Bestion Model for interfacial friction implementation test
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragMultiplierModels/twoPhaseDragMultiplierModel`
      - Author:
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragMultiplierModels/Kaiser88/Kaiser88TwoPhaseDragMultiplier`
      - Model by Kaiser et al
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragMultiplierModels/KottowskiSavatteri/KottowskiSavatteriTwoPhaseDragMultiplier`
      - Model by Kottowski and Savatteri for the multiplier, specifically for
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragMultiplierModels/LockhartMartinelli/LockhartMartinelli`
      - Class to handle fluid-structure friction coefficient (drag) multipliers
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragMultiplierModels/ChenKalish/ChenKalishTwoPhaseDragMultiplier`
      - Model by Chen and Kalish for the multiplier, specifically for
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragMultiplierModels/LottesFlinn/LottesFlinnTwoPhaseDragMultiplier`
      - Model in the spirit of Lottes-Flinn, i
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragMultiplierModels/LottesFlinnNguyen/LottesFlinnNguyenTwoPhaseDragMultiplier`
      - Model in the spirit of Lottes-Flinn, i
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragMultiplierModels/constant/constantTwoPhaseDragMultiplier`
      - Returns a constant multiplier
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/twoPhaseDragMultiplierModels/Kaiser74/Kaiser74TwoPhaseDragMultiplier`
      - Model by Kaiser et al
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/FSDragCoefficientModels/FSDragCoefficientModel`
      - Class to manage the fluid-structure (i
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/FSDragCoefficientModels/Rehme/RehmeFSDragCoefficient`
      - Rheme correlation for pressure drop
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/FSDragCoefficientModels/Churchill/ChurchillFSDragCoefficient`
      - Churcill correlation for pressure drop
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/FSDragCoefficientModels/Colebrook/ColebrookFSDragCoefficient`
      - Drag coefficient in the form
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/FSDragCoefficientModels/NoKazimi/NoKazimiFSDragCoefficient`
      - NoKazimi correlation for fluid-structure pressure drop
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/FSDragCoefficientModels/modifiedEngel/modifiedEngelFSDragCoefficient`
      - modifiedEngel correlation for pressure drop
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/FSDragCoefficientModels/BaxiDalleDonne/BaxiDalleDonneFSDragCoefficient`
      - Baxi-Dalle Donne correlation for pressure drop
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/FSDragCoefficientModels/Engel/EngelFSDragCoefficient`
      - Engel correlation for pressure drop
    * - :doc:`generated/porousMediaModels/physicsModels/dragModels/FSDragCoefficientModels/ReynoldsPower/ReynoldsPowerFSDragCoefficient`
      - Drag coefficient in the form


Heat transfer Models
^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/FSHeatTransferCoefficientModel
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/NusseltWallAndHfromFMU/NusseltWallAndHfromFMUFSHeatTransferCoefficient
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/Gorenflo/GorenfloFSHeatTransferCoefficient
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/multiRegimeBoilingTRACE/multiRegimeBoilingTRACEFSHeatTransferCoefficient
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/multiRegimeBoiling/multiRegimeBoilingFSHeatTransferCoefficient
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/Shah/ShahFSHeatTransferCoefficient
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/Nusselt/NusseltFSHeatTransferCoefficient
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/NusseltAndWall/NusseltAndWallFSHeatTransferCoefficient
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/multiRegimeBoilingVapourTRACE/multiRegimeBoilingVapourTRACEFSHeatTransferCoefficient
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/subCooledBoilingFractionModels/subCooledBoilingFractionModel
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/subCooledBoilingFractionModels/SahaZuber/SahaZuberSubCooledBoilingFraction
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/subCooledBoilingFractionModels/constant/constantSubCooledBoilingFraction
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/postCHFModels/CachardLiquid/CachardLiquidFSHeatTransferCoefficient
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/postCHFModels/CachardVapour/CachardVapourFSHeatTransferCoefficient
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/TONBModels/TONBModel
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/TONBModels/Basu/BasuTONB
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/suppressionFactorModels/suppressionFactorModel
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/suppressionFactorModels/Chen/ChenSuppressionFactor
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/suppressionFactorModels/COBRA-TF/COBRA-TFSuppressionFactor
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/CHFModels/criticalHeatFluxModels/CHFModel
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/CHFModels/criticalHeatFluxModels/lookUpTableCHF/lookUpTableCHF
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/CHFModels/criticalHeatFluxModels/constantCHF/constantCHF
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/CHFModels/LeidenFrostTemperatureModels/TLFModel
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/CHFModels/LeidenFrostTemperatureModels/GroeneveldStewart/GroeneveldStewartTLF
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/superpositionNucleateBoiling/superpositionNucleateBoilingFSHeatTransferCoefficient
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/flowEnhancementFactorModels/flowEnhancementFactorModel
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/flowEnhancementFactorModels/RezkallahSims/RezkallahSimsFlowEnhancementFactor
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/flowEnhancementFactorModels/Chen/ChenFlowEnhancementFactor
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/flowEnhancementFactorModels/COBRA-TF/COBRA-TFFlowEnhancementFactor
   generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/multiRegimeBoilingTRACECHF/multiRegimeBoilingTRACECHFFSHeatTransferCoefficient
   generated/porousMediaModels/physicsModels/heatTransferModels/FFHeatTransferCoefficientModels/FFHeatTransferCoefficientModel
   generated/porousMediaModels/physicsModels/heatTransferModels/FFHeatTransferCoefficientModels/NoKazimi/NoKazimiFFHeatTransferCoefficient
   generated/porousMediaModels/physicsModels/heatTransferModels/FFHeatTransferCoefficientModels/Nusselt/NusseltFFHeatTransferCoefficient

.. list-table::
    :widths: 50 50

    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/FSHeatTransferCoefficientModel`
      - Class to manage the fluid-structure (i
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/NusseltWallAndHfromFMU/NusseltWallAndHfromFMUFSHeatTransferCoefficient`
      - Heat transfer coefficient for forced convective flows that is computed
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/Gorenflo/GorenfloFSHeatTransferCoefficient`
      - Heat transfer coefficient for pool boiling scenarios calculated via the
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/multiRegimeBoilingTRACE/multiRegimeBoilingTRACEFSHeatTransferCoefficient`
      - Heat transfer coefficient for that covers many regimes, from single phase
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/multiRegimeBoiling/multiRegimeBoilingFSHeatTransferCoefficient`
      - Heat transfer coefficient that covers many regimes, from single phase
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/Shah/ShahFSHeatTransferCoefficient`
      - Heat transfer coefficient for pool boiling scenarios calculated via the
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/Nusselt/NusseltFSHeatTransferCoefficient`
      - Heat transfer coefficient for forced convective flows that is computed
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/NusseltAndWall/NusseltAndWallFSHeatTransferCoefficient`
      - Heat transfer coefficient for forced convective flows that is computed
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/multiRegimeBoilingVapourTRACE/multiRegimeBoilingVapourTRACEFSHeatTransferCoefficient`
      - Heat transfer coefficient for the transition Boiling Regime and the
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/subCooledBoilingFractionModels/subCooledBoilingFractionModel`
      - Model for the prediction of the sub-cooled bulk liquid temperature at
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/subCooledBoilingFractionModels/SahaZuber/SahaZuberSubCooledBoilingFraction`
      - Calculate the temperature of onset of nucleate boiling according to the
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/subCooledBoilingFractionModels/constant/constantSubCooledBoilingFraction`
      - Return constant subCooled boiling fraction
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/postCHFModels/CachardLiquid/CachardLiquidFSHeatTransferCoefficient`
      - Heat transfer coefficient for Inverted Annular Flow - for the liquid phase
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/postCHFModels/CachardVapour/CachardVapourFSHeatTransferCoefficient`
      - Heat transfer coefficient for Inverted Annular Flow - for the vapour phase
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/TONBModels/TONBModel`
      - Model for the prediction of the temperature of the onset of nucleate
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/TONBModels/Basu/BasuTONB`
      - Calculate the temperature of onset of nucleate boiling according to the
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/suppressionFactorModels/suppressionFactorModel`
      - Author:
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/suppressionFactorModels/Chen/ChenSuppressionFactor`
      - Author:
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/suppressionFactorModels/COBRA-TF/COBRA-TFSuppressionFactor`
      - Author:
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/CHFModels/criticalHeatFluxModels/CHFModel`
      - Model for the prediction of the critical heat flux
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/CHFModels/criticalHeatFluxModels/lookUpTableCHF/lookUpTableCHF`
      - Compute the Critical Heat Flux with a Look Up Table
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/CHFModels/criticalHeatFluxModels/constantCHF/constantCHF`
      - Gives a constant Critical Heat Flux provided in Phase Properties
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/CHFModels/LeidenFrostTemperatureModels/TLFModel`
      - Model for the prediction of the Leidenfrost temperature (i
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/CHFModels/LeidenFrostTemperatureModels/GroeneveldStewart/GroeneveldStewartTLF`
      - Calculate the Leidenfrost Temperature at which the heat flux is minimal in the
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/superpositionNucleateBoiling/superpositionNucleateBoilingFSHeatTransferCoefficient`
      - Heat transfer coefficient for both single-phase convective scenarios that
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/flowEnhancementFactorModels/flowEnhancementFactorModel`
      - Author:
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/flowEnhancementFactorModels/RezkallahSims/RezkallahSimsFlowEnhancementFactor`
      - Two-phase flow enhancement factor for the convective heat transfer
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/flowEnhancementFactorModels/Chen/ChenFlowEnhancementFactor`
      - Author:
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/subModels/flowEnhancementFactorModels/COBRA-TF/COBRA-TFFlowEnhancementFactor`
      - Author:
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FSHeatTransferCoefficientModels/multiRegimeBoilingTRACECHF/multiRegimeBoilingTRACECHFFSHeatTransferCoefficient`
      - Heat transfer coefficient for that covers many regimes, from single phase
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FFHeatTransferCoefficientModels/FFHeatTransferCoefficientModel`
      - Class to manage the fluid-structure (i
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FFHeatTransferCoefficientModels/NoKazimi/NoKazimiFFHeatTransferCoefficient`
      - NoKazimi model for the liquid side of a liquid-vapour heat transfer
    * - :doc:`generated/porousMediaModels/physicsModels/heatTransferModels/FFHeatTransferCoefficientModels/Nusselt/NusseltFFHeatTransferCoefficient`
      - Heat transfer coefficient for forced convective flows that is computed


Regime Map Models
^^^^^^^^^^^^^^^^^

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/porousMediaModels/physicsModels/regimeMapModels/regimeMapModel
   generated/porousMediaModels/physicsModels/regimeMapModels/twoParameters/regimeBoundary2D
   generated/porousMediaModels/physicsModels/regimeMapModels/twoParameters/twoParameters
   generated/porousMediaModels/physicsModels/regimeMapModels/twoParameters/regimeDomain2D
   generated/porousMediaModels/physicsModels/regimeMapModels/oneParameter/oneParameter

.. list-table::
    :widths: 50 50

    * - :doc:`generated/porousMediaModels/physicsModels/regimeMapModels/regimeMapModel`
      - Run-time selectable class to handle flow regime map models
    * - :doc:`generated/porousMediaModels/physicsModels/regimeMapModels/twoParameters/regimeBoundary2D`
      - 2D line segment representing the boundary between two different
    * - :doc:`generated/porousMediaModels/physicsModels/regimeMapModels/twoParameters/twoParameters`
      - Regime map model that depends on two parameters
    * - :doc:`generated/porousMediaModels/physicsModels/regimeMapModels/twoParameters/regimeDomain2D`
      - Regime domain representation in 2D parameter-space
    * - :doc:`generated/porousMediaModels/physicsModels/regimeMapModels/oneParameter/oneParameter`
      - Regime map model that depend on a single scalar parameter value


Contact Partition Models
^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/porousMediaModels/physicsModels/contactPartitionModels/contactPartitionModel
   generated/porousMediaModels/physicsModels/contactPartitionModels/complementary/complementaryContactPartition
   generated/porousMediaModels/physicsModels/contactPartitionModels/linear/linearContactPartition

.. list-table::
    :widths: 50 50

    * - :doc:`generated/porousMediaModels/physicsModels/contactPartitionModels/contactPartitionModel`
      - Class to handle fluid-structure interfacial contact contact fractions
    * - :doc:`generated/porousMediaModels/physicsModels/contactPartitionModels/complementary/complementaryContactPartition`
      - Author:
    * - :doc:`generated/porousMediaModels/physicsModels/contactPartitionModels/linear/linearContactPartition`
      - Author:


Dispersion Models
^^^^^^^^^^^^^^^^^

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/porousMediaModels/physicsModels/dispersionModels/dispersionModel
   generated/porousMediaModels/physicsModels/dispersionModels/constant/constantDispersion

.. list-table::
    :widths: 50 50

    * - :doc:`generated/porousMediaModels/physicsModels/dispersionModels/dispersionModel`
      - Class to handle fluid-structure interfacial contact contact fractions
    * - :doc:`generated/porousMediaModels/physicsModels/dispersionModels/constant/constantDispersion`
      - Class to handle fluid-structure interfacial contact contact fractions


Fluid Diameter Models
^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/porousMediaModels/physicsModels/fluidDiameterModels/fluidDiameterModel
   generated/porousMediaModels/physicsModels/fluidDiameterModels/isomolarBubble/isomolarBubbleFluidDiameter
   generated/porousMediaModels/physicsModels/fluidDiameterModels/pipeFilm/pipeFilmFluidDiameter
   generated/porousMediaModels/physicsModels/fluidDiameterModels/isothermalBubble/isothermalBubbleFluidDiameter
   generated/porousMediaModels/physicsModels/fluidDiameterModels/WallisFilm/WallisFilmFluidDiameter

.. list-table::
    :widths: 50 50

    * - :doc:`generated/porousMediaModels/physicsModels/fluidDiameterModels/fluidDiameterModel`
      - Class to handle fluid diameter
    * - :doc:`generated/porousMediaModels/physicsModels/fluidDiameterModels/isomolarBubble/isomolarBubbleFluidDiameter`
      - Compute bubble/droplet diameter based on pressure and temperature
    * - :doc:`generated/porousMediaModels/physicsModels/fluidDiameterModels/pipeFilm/pipeFilmFluidDiameter`
      - Correlation to be used in annular flow scenarios in which the phase of
    * - :doc:`generated/porousMediaModels/physicsModels/fluidDiameterModels/isothermalBubble/isothermalBubbleFluidDiameter`
      - Compute bubble/droplet diameter based on pressure changes assuming
    * - :doc:`generated/porousMediaModels/physicsModels/fluidDiameterModels/WallisFilm/WallisFilmFluidDiameter`
      - Correlation to be used in annular flow scenarios in which the phase of


Interfacial Area Models
^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/porousMediaModels/physicsModels/interfacialAreaModels/interfacialAreaModel
   generated/porousMediaModels/physicsModels/interfacialAreaModels/NoKazimi/NoKazimiInterfacialArea
   generated/porousMediaModels/physicsModels/interfacialAreaModels/annular/annularInterfacialArea
   generated/porousMediaModels/physicsModels/interfacialAreaModels/Schor/SchorInterfacialArea
   generated/porousMediaModels/physicsModels/interfacialAreaModels/spherical/sphericalInterfacialArea

.. list-table::
    :widths: 50 50

    * - :doc:`generated/porousMediaModels/physicsModels/interfacialAreaModels/interfacialAreaModel`
      - Class to handle fluid-fluid interfacial area
    * - :doc:`generated/porousMediaModels/physicsModels/interfacialAreaModels/NoKazimi/NoKazimiInterfacialArea`
      - Author:
    * - :doc:`generated/porousMediaModels/physicsModels/interfacialAreaModels/annular/annularInterfacialArea`
      - Interfacial area model that assumes that the dispersed phase flow is
    * - :doc:`generated/porousMediaModels/physicsModels/interfacialAreaModels/Schor/SchorInterfacialArea`
      - Author:
    * - :doc:`generated/porousMediaModels/physicsModels/interfacialAreaModels/spherical/sphericalInterfacialArea`
      - Author:


Phase Change Models
^^^^^^^^^^^^^^^^^^^

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/porousMediaModels/physicsModels/phaseChangeModels/phaseChangeModel
   generated/porousMediaModels/physicsModels/phaseChangeModels/latentHeatModels/latentHeatModel
   generated/porousMediaModels/physicsModels/phaseChangeModels/latentHeatModels/water/waterLatentHeat
   generated/porousMediaModels/physicsModels/phaseChangeModels/latentHeatModels/FinkLeibowitz/FinkLeibowitzLatentHeat
   generated/porousMediaModels/physicsModels/phaseChangeModels/latentHeatModels/fromThermophysicalProperties/fromThermophysicalPropertiesLatentHeat
   generated/porousMediaModels/physicsModels/phaseChangeModels/heatDriven/heatDrivenPhaseChange
   generated/porousMediaModels/physicsModels/phaseChangeModels/saturationModels/saturationModel
   generated/porousMediaModels/physicsModels/phaseChangeModels/saturationModels/water/waterSaturation
   generated/porousMediaModels/physicsModels/phaseChangeModels/saturationModels/BrowningPotter/BrowningPotterSaturation
   generated/porousMediaModels/physicsModels/phaseChangeModels/saturationModels/constantTemperature/constantTemperatureSaturation
   generated/porousMediaModels/physicsModels/phaseChangeModels/saturationModels/waterTRACE/waterTRACESaturation
   generated/porousMediaModels/physicsModels/phaseChangeModels/forcedConstant/forcedConstantPhaseChange

.. list-table::
    :widths: 50 50

    * - :doc:`generated/porousMediaModels/physicsModels/phaseChangeModels/phaseChangeModel`
      - SourceFiles
    * - :doc:`generated/porousMediaModels/physicsModels/phaseChangeModels/latentHeatModels/latentHeatModel`
      - Class used to describe the latent heat of fluid that is undergoing phase
    * - :doc:`generated/porousMediaModels/physicsModels/phaseChangeModels/latentHeatModels/water/waterLatentHeat`
      - This class describes the latent heat of vaporization of water in the
    * - :doc:`generated/porousMediaModels/physicsModels/phaseChangeModels/latentHeatModels/FinkLeibowitz/FinkLeibowitzLatentHeat`
      - This class describes the latent heat of vaporization of liquid sodium
    * - :doc:`generated/porousMediaModels/physicsModels/phaseChangeModels/latentHeatModels/fromThermophysicalProperties/fromThermophysicalPropertiesLatentHeat`
      - Class that computes latent heat based on the specified enthalpies of
    * - :doc:`generated/porousMediaModels/physicsModels/phaseChangeModels/heatDriven/heatDrivenPhaseChange`
      - Phase change models that compute mass transfer based on interfacial
    * - :doc:`generated/porousMediaModels/physicsModels/phaseChangeModels/saturationModels/saturationModel`
      - SourceFiles
    * - :doc:`generated/porousMediaModels/physicsModels/phaseChangeModels/saturationModels/water/waterSaturation`
      - Saturation model for water in the 0
    * - :doc:`generated/porousMediaModels/physicsModels/phaseChangeModels/saturationModels/BrowningPotter/BrowningPotterSaturation`
      - Saturation model based on results by Browning and Potter, refer to
    * - :doc:`generated/porousMediaModels/physicsModels/phaseChangeModels/saturationModels/constantTemperature/constantTemperatureSaturation`
      - Constant saturation pressure and temperature
    * - :doc:`generated/porousMediaModels/physicsModels/phaseChangeModels/saturationModels/waterTRACE/waterTRACESaturation`
      - Saturation model for water in the 0
    * - :doc:`generated/porousMediaModels/physicsModels/phaseChangeModels/forcedConstant/forcedConstantPhaseChange`
      - Phase change models that compute mass transfer based on interfacial


Phase Pairs
^^^^^^^^^^^

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/porousMediaModels/physicsModels/phasePairs/FFPair
   generated/porousMediaModels/physicsModels/phasePairs/FSPair

.. list-table::
    :widths: 50 50

    * - :doc:`generated/porousMediaModels/physicsModels/phasePairs/FFPair`
      - Class to handle dimensionless numbers that characterise the system, e
    * - :doc:`generated/porousMediaModels/physicsModels/phasePairs/FSPair`
      - Class to handle dimensionless numbers that characterise the system, e


Templated Models
^^^^^^^^^^^^^^^^

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/porousMediaModels/physicsModels/templatedModels/byRegime/byRegimeModel
   generated/porousMediaModels/physicsModels/templatedModels/constant/constantModel

.. list-table::
    :widths: 50 50

    * - :doc:`generated/porousMediaModels/physicsModels/templatedModels/byRegime/byRegimeModel`
      - Author:
    * - :doc:`generated/porousMediaModels/physicsModels/templatedModels/constant/constantModel`
      - Author:


Turbulence Models
^^^^^^^^^^^^^^^^^

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/porousMediaModels/physicsModels/turbulenceModels/mixtureKEpsilon/mixtureKEpsilon
   generated/porousMediaModels/physicsModels/turbulenceModels/porousKEpsilon2PhaseCorrected/porousKEpsilon2PhaseCorrected
   generated/porousMediaModels/physicsModels/turbulenceModels/porousKEpsilon/porousKEpsilon
   generated/porousMediaModels/physicsModels/turbulenceModels/LaheyKEpsilon/LaheyKEpsilon

.. list-table::
    :widths: 50 50

    * - :doc:`generated/porousMediaModels/physicsModels/turbulenceModels/mixtureKEpsilon/mixtureKEpsilon`
      - Mixture k-epsilon turbulence model for two-phase gas-liquid systems
    * - :doc:`generated/porousMediaModels/physicsModels/turbulenceModels/porousKEpsilon2PhaseCorrected/porousKEpsilon2PhaseCorrected`
      - Same as porousKEpsiolon but turbulent intensity correlation in the form
    * - :doc:`generated/porousMediaModels/physicsModels/turbulenceModels/porousKEpsilon/porousKEpsilon`
      - Same as standard OpenFOAM
    * - :doc:`generated/porousMediaModels/physicsModels/turbulenceModels/LaheyKEpsilon/LaheyKEpsilon`
      - Continuous-phase k-epsilon model including bubble-generated turbulence


Virutal Mass Models
^^^^^^^^^^^^^^^^^^^

.. toctree::
   :hidden:
   :maxdepth: 1


.. list-table::
    :widths: 50 50



-----------------------------------------
Structural mechanics and fuel performance
-----------------------------------------



Physics Sub-Solvers
-------------------

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/offbeat/src/offbeatLib/physicsSubSolvers/flowSubSolver/rhoPimpleFlowSubSolver
   generated/offbeat/src/offbeatLib/physicsSubSolvers/flowSubSolver/pimpleFlowSubSolver
   generated/offbeat/src/offbeatLib/physicsSubSolvers/flowSubSolver/flowSubSolver
   generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/largeStrainUpdLag
   generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/smallStrainIncUpdated
   generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/largeStrainTotLag
   generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/smallStrain
   generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/mechanicsSubSolver
   generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/mechanicsSubSolverTemplates
   generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/multiMaterialInterface/topoSetSourceMultiMaterialInterface
   generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/multiMaterialInterface/uniformDirectionalMultiMaterialInterface
   generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/multiMaterialInterface/cellZoneMultiMaterialInterface
   generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/multiMaterialInterface/uniformMultiMaterialInterface
   generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/multiMaterialInterface/uniform2DMultiMaterialInterface
   generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/multiMaterialInterface/multiMaterialInterface
   generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/multiMaterialInterface/faceSetMultiMaterialInterface
   generated/offbeat/src/offbeatLib/physicsSubSolvers/neutronicsSubSolver/neutronicsSubSolver
   generated/offbeat/src/offbeatLib/physicsSubSolvers/neutronicsSubSolver/diffusionSolver
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/elementTransportByList
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/elementTransport
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/PuRedistributionSolver
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/porosityTransportSolver
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/AmRedistributionSolver
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/fissionProductsDiffusionSolver
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/porosityTransportDynMeshSolver
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/transportSolver
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HSolverOneTSS
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HNGD
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HSolver
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HSolverEq
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/hydrogenTransportSolver
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSdModel/TSSdDXD4Liner
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSdModel/TSSdModel
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSdModel/TSSdBISON
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSdModel/TSSdUneZrLiner
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSdModel/TSSdUneZr2
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSdModel/TSSdMcMinn
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSdModel/TSSdKammenzind
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSdModel/TSSdGeneric
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSdModel/TSSdMerlino
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/hydrideReorientationModel/hydrideReorientationModel
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/hydrideReorientationModel/hydrideReorientationDesquines
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HDiffusionCoefficientModel/HDiffusionCoefficientMallet
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HDiffusionCoefficientModel/HDiffusionCoefficientModel
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HDiffusionCoefficientModel/HDiffusionCoefficientKearns
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HDiffusionCoefficientModel/HDiffusionCoefficientGeneric
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HDiffusionCoefficientModel/HDiffusionCoefficientKammenzind
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HDiffusionCoefficientModel/HDiffusionCoefficientBISON
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpGeneric
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpDXD4Liner
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpKammenzind
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpBISON
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpUneZrLiner
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpMerlinoValue
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpModel
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpUneZr2
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpMcMinn
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpMerlinoExp
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/poreVelocityModel/poreVelocityConstant
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/poreVelocityModel/poreVelocityMOXLackey
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/poreVelocityModel/poreVelocityModel
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/poreVelocityModel/poreVelocityMOXVaporPressure
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/poreVelocityModel/poreVelocityUO2Sens
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/poreVelocityModel/poreVelocityMOXClementFinnis
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/diffusionCoefficient/diffCoefPyC
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/diffusionCoefficient/diffCoefKernel
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/diffusionCoefficient/diffCoefModel
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/diffusionCoefficient/diffCoefSiC
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/diffusionCoefficient/diffCoefArrhenius
   generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/diffusionCoefficient/diffCoefGraphite
   generated/offbeat/src/offbeatLib/physicsSubSolvers/thermalSubSolver/solidConductionSolver
   generated/offbeat/src/offbeatLib/physicsSubSolvers/thermalSubSolver/readTemperature
   generated/offbeat/src/offbeatLib/physicsSubSolvers/thermalSubSolver/thermalSubSolver

.. list-table::
    :widths: 50 50

    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/flowSubSolver/rhoPimpleFlowSubSolver`
      - Flow solver based using PIMPLE algorithm
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/flowSubSolver/pimpleFlowSubSolver`
      - Flow solver based using PIMPLE algorithm
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/flowSubSolver/flowSubSolver`
      - Base class for flow solvers
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/largeStrainUpdLag`
      - This class olves for incremental displacement `DD` in an updated Lagrangian configuration 
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/smallStrainIncUpdated`
      - This class olves for incremental displacement `DD` in an updated configuration 
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/largeStrainTotLag`
      - This class olves for total displacement `D` in a total Lagrangian configuration 
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/smallStrain`
      - This class solves for total displacement `D` in a total Lagrangian configuration 
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/mechanicsSubSolver`
      - When using the `constant` mechanics solver class, the displacement field `D`
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/mechanicsSubSolverTemplates`
      - No description available
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/multiMaterialInterface/topoSetSourceMultiMaterialInterface`
      - This class applies a multi-material correction to the solution of the 
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/multiMaterialInterface/uniformDirectionalMultiMaterialInterface`
      - Interfaces between materials are defined by a directionally-
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/multiMaterialInterface/cellZoneMultiMaterialInterface`
      - This class applies a multi-material correction to the solution of the 
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/multiMaterialInterface/uniformMultiMaterialInterface`
      - This class applies a uniform multi-material correction to the solution of the 
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/multiMaterialInterface/uniform2DMultiMaterialInterface`
      - This class applies a uniform multi-material correction to the solution of the 
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/multiMaterialInterface/multiMaterialInterface`
      - Base class for defining the set of faces that constitute the interfaces
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/mechanicsSubSolver/multiMaterialInterface/faceSetMultiMaterialInterface`
      - This class applies a multi-material correction to the solution of the 
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/neutronicsSubSolver/neutronicsSubSolver`
      - As the name implies, the `none` neutronics solver class in OFFBEAT 
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/neutronicsSubSolver/diffusionSolver`
      - In the `diffusionSolver` solver, the neutron flux distribution is 
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/elementTransportByList`
      - Class that gathers all the user-selected transport solvers in a ptrList, 
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/elementTransport`
      - Mother class for handling element transport solvers
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/PuRedistributionSolver`
      - Daugther class of transportSolver class that solves for the redistribution
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/porosityTransportSolver`
      - Daugther class of transportSolver class that solves for the redistribution
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/AmRedistributionSolver`
      - Daugther class of transportSolver class that solves for the redistribution
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/fissionProductsDiffusionSolver`
      - Solver for Fission products transport
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/porosityTransportDynMeshSolver`
      - Daugther class of transportSolver class that solves for the redistribution
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/transportSolver`
      - Mother class for transport solvers
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HSolverOneTSS`
      - No description available
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HNGD`
      - No description available
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HSolver`
      - No description available
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HSolverEq`
      - No description available
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/hydrogenTransportSolver`
      - Mother class for the hydrogen transport solvers
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSdModel/TSSdDXD4Liner`
      - TSSD model for the DXD4 outer liner
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSdModel/TSSdModel`
      - yield stress mother class
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSdModel/TSSdBISON`
      - varying yield stress as a function of plastic strain
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSdModel/TSSdUneZrLiner`
      - TSSD model for the pure Zirconium inner liner
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSdModel/TSSdUneZr2`
      - TSSD model for the Zr-2 substrate (used in Inner liner claddings)
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSdModel/TSSdMcMinn`
      - varying yield stress as a function of plastic strain
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSdModel/TSSdKammenzind`
      - varying yield stress as a function of plastic strain
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSdModel/TSSdGeneric`
      - varying yield stress as a function of plastic strain
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSdModel/TSSdMerlino`
      - varying yield stress as a function of plastic strain
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/hydrideReorientationModel/hydrideReorientationModel`
      - No description available
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/hydrideReorientationModel/hydrideReorientationDesquines`
      - No description available
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HDiffusionCoefficientModel/HDiffusionCoefficientMallet`
      - varying yield stress as a function of plastic strain
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HDiffusionCoefficientModel/HDiffusionCoefficientModel`
      - Base class for the diffusion coefficient of hydrogen ions
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HDiffusionCoefficientModel/HDiffusionCoefficientKearns`
      - varying yield stress as a function of plastic strain
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HDiffusionCoefficientModel/HDiffusionCoefficientGeneric`
      - Usage
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HDiffusionCoefficientModel/HDiffusionCoefficientKammenzind`
      - Diffusion coefficient of Hydrogen in Zircaloy-4 from Kammenzind et al
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/HDiffusionCoefficientModel/HDiffusionCoefficientBISON`
      - varying yield stress as a function of plastic strain
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpGeneric`
      - varying yield stress as a function of plastic strain
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpDXD4Liner`
      - TSSP model for the DXD4 outer liner
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpKammenzind`
      - varying yield stress as a function of plastic strain
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpBISON`
      - varying yield stress as a function of plastic strain
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpUneZrLiner`
      - TSSP model for the pure Zirconium inner liner
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpMerlinoValue`
      - varying yield stress as a function of plastic strain
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpModel`
      - yield stress mother class
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpUneZr2`
      - TSSP model for the Zr-2 substrate (used in Inner liner claddings)
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpMcMinn`
      - varying yield stress as a function of plastic strain
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/hydrogenTransport/TSSpModel/TSSpMerlinoExp`
      - varying yield stress as a function of plastic strain
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/poreVelocityModel/poreVelocityConstant`
      - By selecting this constant model the velocity is read from the previous time
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/poreVelocityModel/poreVelocityMOXLackey`
      - Model to compute pore velocity within MOX fuel with the correlation provided
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/poreVelocityModel/poreVelocityModel`
      - Mother class for pore velocity models
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/poreVelocityModel/poreVelocityMOXVaporPressure`
      - Usage
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/poreVelocityModel/poreVelocityUO2Sens`
      - Model to compute pore velocity with the correlation provided by Sens (P
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/transportSolvers/poreVelocityModel/poreVelocityMOXClementFinnis`
      - Model to compute pore velocity with the correlation provided by Clement and
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/diffusionCoefficient/diffCoefPyC`
      - Diffusion coefficient for PyC material from PARFUME manual and INL report
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/diffusionCoefficient/diffCoefKernel`
      - Diffusion coefficient for Kernel material from PARFUME manual and INL report
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/diffusionCoefficient/diffCoefModel`
      - Mother class for fission products diffusion coefficient models
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/diffusionCoefficient/diffCoefSiC`
      - Diffusion coefficient for PyC material from PARFUME manual and INL report
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/diffusionCoefficient/diffCoefArrhenius`
      - Diffusion coefficient to be calculated in Arrhenius form
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/elementTransport/diffusionCoefficient/diffCoefGraphite`
      - Diffusion coefficient for PyC material from PARFUME manual and INL report
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/thermalSubSolver/solidConductionSolver`
      - In the `solidConduction` solver, the temperature distribution is obtained from
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/thermalSubSolver/readTemperature`
      - The `readTemperature` thermal solver class in OFFBEAT reads the temperature `T`
    * - :doc:`generated/offbeat/src/offbeatLib/physicsSubSolvers/thermalSubSolver/thermalSubSolver`
      - When using the `constant` thermal solver class, the temperature field `T`


Heat Sources
------------

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/offbeat/src/offbeatLib/heatSource/heatSource
   generated/offbeat/src/offbeatLib/heatSource/constantLhgr
   generated/offbeat/src/offbeatLib/heatSource/timeDependentVhgr
   generated/offbeat/src/offbeatLib/heatSource/constantHeatSource
   generated/offbeat/src/offbeatLib/heatSource/laserHeatSource
   generated/offbeat/src/offbeatLib/heatSource/fmiLhgr
   generated/offbeat/src/offbeatLib/heatSource/timeDependentLhgr

.. list-table::
    :widths: 50 50

    * - :doc:`generated/offbeat/src/offbeatLib/heatSource/heatSource`
      - As the name implies, the `none` heat source class in OFFBEAT allows you to switch
    * - :doc:`generated/offbeat/src/offbeatLib/heatSource/constantLhgr`
      - The `constantLhgr` heatSource class in OFFBEAT allows you to set a constant
    * - :doc:`generated/offbeat/src/offbeatLib/heatSource/timeDependentVhgr`
      - The `timeDependentVhgr` class in OFFBEAT allows you to set a average
    * - :doc:`generated/offbeat/src/offbeatLib/heatSource/constantHeatSource`
      - The `constant` heatSource class in OFFBEAT allows the user to set the volumetric
    * - :doc:`generated/offbeat/src/offbeatLib/heatSource/laserHeatSource`
      - This heat source model is used to simulate the heating due to energy 
    * - :doc:`generated/offbeat/src/offbeatLib/heatSource/fmiLhgr`
      - Derived from the base heatSource class, this class act as the parent class 
    * - :doc:`generated/offbeat/src/offbeatLib/heatSource/timeDependentLhgr`
      - Similar to the `constantLhgr` heatSource class, the `timeDependentLhgr` class in


Fast Flux
---------

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/offbeat/src/offbeatLib/fastFlux/timeDependentAxialProfile
   generated/offbeat/src/offbeatLib/fastFlux/fastFlux
   generated/offbeat/src/offbeatLib/fastFlux/constantFastFlux

.. list-table::
    :widths: 50 50

    * - :doc:`generated/offbeat/src/offbeatLib/fastFlux/timeDependentAxialProfile`
      - The `timeDependentAxialProfile` class in OFFBEAT allows you to set a time
    * - :doc:`generated/offbeat/src/offbeatLib/fastFlux/fastFlux`
      - As the name implies, the `none` fast flux class in OFFBEAT allows you to switch
    * - :doc:`generated/offbeat/src/offbeatLib/fastFlux/constantFastFlux`
      - The `constant` fastFlux class in OFFBEAT allows you to set the fast flux


Materials
---------

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/offbeat/src/offbeatLib/materials/materialModel/buffer
   generated/offbeat/src/offbeatLib/materials/materialModel/UPuO2
   generated/offbeat/src/offbeatLib/materials/materialModel/materialModel
   generated/offbeat/src/offbeatLib/materials/materialModel/HastelloyN
   generated/offbeat/src/offbeatLib/materials/materialModel/UO2
   generated/offbeat/src/offbeatLib/materials/materialModel/SiC
   generated/offbeat/src/offbeatLib/materials/materialModel/fuelMaterial
   generated/offbeat/src/offbeatLib/materials/materialModel/PyC
   generated/offbeat/src/offbeatLib/materials/materialModel/constantMaterial
   generated/offbeat/src/offbeatLib/materials/materialModel/steel1515Ti
   generated/offbeat/src/offbeatLib/materials/materialModel/molybdenum
   generated/offbeat/src/offbeatLib/materials/materialModel/inconel600
   generated/offbeat/src/offbeatLib/materials/materialModel/zircaloy
   generated/offbeat/src/offbeatLib/materials/materialModel/damageModel/isotropicCracking
   generated/offbeat/src/offbeatLib/materials/materialModel/damageModel/mazarsDamageModel
   generated/offbeat/src/offbeatLib/materials/materialModel/damageModel/damageModel
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/relocation/relocationFRAPCON
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/relocation/relocationModel
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/densification/densificationEmpirical
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/densification/densificationModel
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/densification/densificationFRAPCON
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingPARFUMEPyCdata
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingGrowthMatproZy
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingGrowthGeneralized1515Ti
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingWrightShamHastelloyN
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingGrowthAIM11515Ti
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingMATPRO
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingFRAPCON
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingGrowthBISONZy
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingCorrelationPyC
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/constantSwelling
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingPARFUMEBuffer
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingFrCrAl
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingModel
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingFBRMOX
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingPARFUMEPyC
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/failureModels/ZyRIAJernkvistModified
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/failureModels/Zy_CSED_RIA_EPRI
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/failureModels/UPuO2meltingMagni2020
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/failureModels/failureModel
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/failureModels/UO2meltingMATPRO
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/failureModels/ZyOverstressBISON
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/failureModels/ZyOverstrainBISON
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/failureModels/combinedFailureCriteria
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/failureModels/ZyPlasticInstabilityBISON
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/phaseTransition/phaseTransitionModel
   generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/phaseTransition/phaseTransitionZyDynamic
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityConstant
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityModel
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityMatproZy
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityMatproUPuO2
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacitySneadSiC
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityIAEAZy
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityFinkUPuO2
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityMolybdenum
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityBanerjee1515Ti
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityMatproUO2
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionConstant
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionRelapUO2
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionModel
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionMatproZy
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionGehr1515Ti
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionPARFUMEBuffer
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionSneadSiC
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionMatproUPuO2
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionMartinUPuO2
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionSwindemanHastelloyN
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionMolybdenum
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionPARFUMEPyC
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionMAMOX
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionPARFUMESiC
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionLemehovUPuO2
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityUPuO2LanningBeyer
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityUPuO2Kato
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityUPuO2Brancheria
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityMolybdenum
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityPARFUMEBuffer
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityUO2IFA601
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityUPuO2MagniMAMOX
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityUPuO2Philipponneau
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityPARFUMESiC
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityRelapZy
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityModel
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityTobbe1515Ti
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivitySwindemanHastelloyN
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityNfirUO2
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityMatproUO2
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityConstant
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusPARFUMEPyC
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusMatproZy
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusPARFUMEBuffer
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusMatproUPuO2
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusSckCenUPuO2
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusMatproUO2
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusPARFUMESiC
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusWatrousHastelloyN
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusTobbe1515Ti
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusModel
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusMolybdenum
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusSneadSiC
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusConstant
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusHofmanD9
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/PoissonRatioModel
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/PoissonRatioTobbe1515Ti
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/PoissonRatioConstant
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/constantPoissonRatioUPuO2
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/PoissonRatioMatproZy
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/constantPoissonRatioUO2
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/constantPoissonRatioMolybdenum
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/constantPoissonRatioZy
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/emissivity/emissivityRelapUO2
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/emissivity/constantEmissivityMolybdenum
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/emissivity/emissivityConstant
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/emissivity/emissivityModel
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/emissivity/constantEmissivityZy
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/density/constantDensityUO2
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/density/constantDensityMolybdenum
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/density/densityModel
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/density/densityConstant
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/density/densitySchumann1515Ti
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/density/constantDensityUPuO2
   generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/density/densityIAEAZy

.. list-table::
    :widths: 50 50

    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/buffer`
      - buffer properties taken from PARFUME;
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/UPuO2`
      - UPuO2 material class
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/materialModel`
      - Parent class for all materials
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/HastelloyN`
      - HastelloyN thermomechanical properties are taken from different sources
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/UO2`
      - UO2 material class
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/SiC`
      - SiC properties taken mostly from PARFUME;
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/fuelMaterial`
      - Parent class for all fuel materials (UO2, UPuO2 etc
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/PyC`
      - PyC properties taken from PARFUME;
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/constantMaterial`
      - Class for the user manual input of the thermomechanical properties of the 
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/steel1515Ti`
      - Class modelling the austenitic steel 15-15 Ti cladding material
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/molybdenum`
      - Properties for molybdenum, taken mostly from BISON manual
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/inconel600`
      - Inconel600 thermoMechanical properties, taken mostly from the web
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/zircaloy`
      - Zircaloy properties taken mostly from MATPRO, BISON manual and IAEA nuclear
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/damageModel/isotropicCracking`
      - Usage
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/damageModel/mazarsDamageModel`
      - Damage model for UO2
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/damageModel/damageModel`
      - yield stress mother class
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/relocation/relocationFRAPCON`
      - Class modelling the relocation phoenomenon derived from FRAPCON
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/relocation/relocationModel`
      - Parent class for relocation model
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/densification/densificationEmpirical`
      - Empirical model for densification
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/densification/densificationModel`
      - Parent class for densification model
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/densification/densificationFRAPCON`
      - Class modelling the densification phenomenon derived from FRAPCON
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingPARFUMEPyCdata`
      - No description available
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingGrowthMatproZy`
      - Class modelling the irradiation growth phenomenon according to MATPRO 
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingGrowthGeneralized1515Ti`
      - Class modelling the void swelling growth phenomenon for 15-15 Ti cladding 
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingWrightShamHastelloyN`
      - Class modelling the void swelling phenomenon for HastelloyN cladding
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingGrowthAIM11515Ti`
      - Class modelling the void swelling growth phenomenon for 15-15 Ti cladding 
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingMATPRO`
      - Class modelling swelling phenomenon derived from FRAPCON
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingFRAPCON`
      - Class modelling swelling phenomenon derived from FRAPCON
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingGrowthBISONZy`
      - Class modelling the irradiation growth phenomenon derived from Moose documentation
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingCorrelationPyC`
      - Class handling irradiation-induced dimensional change eigenstrain phenomenon
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/constantSwelling`
      - Class caculating the swelling with a provided constant swelling rate
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingPARFUMEBuffer`
      - Class handling irradiation-induced dimensional change phenomenon for Buffer
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingFrCrAl`
      - Class handling swelling phenomenon for FrCrAl derived from BISON manual
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingModel`
      - Parent class for swelling model
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingFBRMOX`
      - Class modelling swelling strain for FBR MOX fuel from Dienst at al
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/swelling/swellingPARFUMEPyC`
      - Class handling irradiation-induced dimensional change eigenstrain phenomenon
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/failureModels/ZyRIAJernkvistModified`
      - Overstrain criterion to asses clad failure during the PCMI phase of a RIA transient
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/failureModels/Zy_CSED_RIA_EPRI`
      - Clad failure criterion during the PCMI phase of a RIA transient based on the strain energy density (SED)
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/failureModels/UPuO2meltingMagni2020`
      - Model to check failure of UPuO2 based on a limit melting temperature from
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/failureModels/failureModel`
      - Mother class for materials failure criteria
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/failureModels/UO2meltingMATPRO`
      - Model to check failure of UO2 based on a limit melting temperature from
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/failureModels/ZyOverstressBISON`
      - Burst overstress criterion for Zircaloy-4 derived from BISON code
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/failureModels/ZyOverstrainBISON`
      - Burst overstrain criterion for Zircaloy-4 derived from BISON code
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/failureModels/combinedFailureCriteria`
      - Criterion used to combine multiple failure criteria
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/failureModels/ZyPlasticInstabilityBISON`
      - Plastic instability criterion for Zircaloy-4 derived from BISON code
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/phaseTransition/phaseTransitionModel`
      - Parent class for phase transition model
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/behavioralModels/phaseTransition/phaseTransitionZyDynamic`
      - Class for the beta phase transition of Zy-based claddings, dynamic approach
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityConstant`
      - Model for constant heat capacity
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityModel`
      - Mother class for heat capacity models
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityMatproZy`
      - Correlation for Zircaloy heat capacity from Matpro
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityMatproUPuO2`
      - Heat Capacity model for (U,Pu)O2 derived from MATPROv11
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacitySneadSiC`
      - Correlation for SiC heat capacity from Snead et al
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityIAEAZy`
      - Correlation for Zircaloy heat capacity from IAEA
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityFinkUPuO2`
      - Heat Capacity model for (U,Pu)O2 derived from Fink
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityMolybdenum`
      - SourceFiles
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityBanerjee1515Ti`
      - Correlation for 15-15 Ti heat capacity from Banerjee et al
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/heatCapacityMatproUO2`
      - Heat Capacity model for UO2 derived from MATPROv11
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionConstant`
      - Model for constant thermal expansion coefficient
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionRelapUO2`
      - Class modelling thermal expansion of UO2 fuel derived from RELAP
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionModel`
      - Mother class for thermal expansion models
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionMatproZy`
      - Class modelling thermal expansion of UO2 fuel from MATPROv11
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionGehr1515Ti`
      - Class modelling the thermal expansion of 15-15 Ti cladding through the
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionPARFUMEBuffer`
      - Class modelling the thermal expansion of Buffer through the code PARFUME
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionSneadSiC`
      - Class modelling the thermal expansion of SiC with a value from Snead
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionMatproUPuO2`
      - Class modelling thermal expansion of (U,Pu)O2 MOX fuel from MATPROv11
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionMartinUPuO2`
      - Class modelling thermal expansion of MOX fuel according to Martin 
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionSwindemanHastelloyN`
      - Class modelling thermal expansion of H-N from Swinderman
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionMolybdenum`
      - Class modelling thermal expansion of Molybdenum
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionPARFUMEPyC`
      - Class modelling the transversly isotropic thermal expansion of PyC through
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionMAMOX`
      - Class modelling  isotropic thermal expansion in MA-MOX fuel (for Pu = 0
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionPARFUMESiC`
      - Class modelling the thermal expansion of SiC with a constant value 4
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/thermalExpansionLemehovUPuO2`
      - Class modelling thermal expansion of (U,Pu)O2 MOX fuel from Lemehov(2020)
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityUPuO2LanningBeyer`
      - SourceFiles
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityUPuO2Kato`
      - S
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityUPuO2Brancheria`
      - Conductivity correlation for MOX fuel Pu=0
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityMolybdenum`
      - Class for Molybdenum conductivity
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityPARFUMEBuffer`
      - Buffer thermal conductivity model
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityUO2IFA601`
      - Model for conductivity of UO2 - IFA-601 rod
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityUPuO2MagniMAMOX`
      - Model for thermal conductivity of high-Pu content MOX
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityUPuO2Philipponneau`
      - Model for thermal conductivity of MOX taken from "Thermal conductivity of 
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityPARFUMESiC`
      - Conductivity of PyC material from PARFUME manual
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityRelapZy`
      - Model for conductivity of Zircaloy material derived by RELAP
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityModel`
      - Mother class for conductivity models
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityTobbe1515Ti`
      - Class modelling conductivity of 15-15 Ti cladding material through the 
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivitySwindemanHastelloyN`
      - Model for conductivity of H-N, derived from the work of Swinderman
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityNfirUO2`
      - Model for thermal conductivity of UO2 from Nfir
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityMatproUO2`
      - Model for UO2 conductivity from MATPROv11
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/conductivityConstant`
      - Model for constant conductivity
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusPARFUMEPyC`
      - Class modelling Young Modulus of PyC from PARFUME
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusMatproZy`
      - Class modelling Young Modulus of Zy from MATPROv11
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusPARFUMEBuffer`
      - Class modelling Young Modulus of Bufffer from PARFUME
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusMatproUPuO2`
      - Class modelling the Young's modulus of (U,Pu)O2 MOX fuel
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusSckCenUPuO2`
      - Correlation for Young's modulus of (U,Pu)O2 MOX fuel from SCK-CEN
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusMatproUO2`
      - Class modelling Young Modulus of UO2 fuel from MATPROv11
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusPARFUMESiC`
      - Class modelling Young Modulus of SiC with the Interpolation from PARFUME
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusWatrousHastelloyN`
      - Class modelling Young Modulus of HN from work of Watrous
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusTobbe1515Ti`
      - Class modelling the Young's modulus of 15-15 Ti using the Tobbe correlation 
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusModel`
      - Mother class for Young's Modulus models
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusMolybdenum`
      - Class modelling Young Modulus of Molybdenum from Bison manual
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusSneadSiC`
      - Class modelling Young Modulus of SiC from Snead et al
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusConstant`
      - Model for constant Young modulus
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/YoungModulusHofmanD9`
      - Class modelling the Young's modulus of D9 steel cladding using the Hofman 
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/PoissonRatioModel`
      - Mother's class for Poisson's Ratio models
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/PoissonRatioTobbe1515Ti`
      - Class modelling Poisson's ratio of 15-15 Ti cladding material based on Tobbe
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/PoissonRatioConstant`
      - Model for constant Poisson Ratio
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/constantPoissonRatioUPuO2`
      - Class modelling the constant Poisson's ratio of UPuO2 MOX fuel from 
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/PoissonRatioMatproZy`
      - Model for Zy Poisson's Ratio from MATPROv11
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/constantPoissonRatioUO2`
      - Class to set UO2 Poisson's Ratio to a constant input value
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/constantPoissonRatioMolybdenum`
      - Class to set Mo Poisson's Ratio to a constant input value
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/constantPoissonRatioZy`
      - Class to set Zy Poisson's Ratio to a constant input value
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/emissivity/emissivityRelapUO2`
      - UO2 emissivity model derived from RELAP
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/emissivity/constantEmissivityMolybdenum`
      - Class to set Mo emissivity to a constant input value
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/emissivity/emissivityConstant`
      - Model for constant emissivity
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/emissivity/emissivityModel`
      - Mother class for emissivity models
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/emissivity/constantEmissivityZy`
      - Class to set Zy emissivity to a constant input value
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/density/constantDensityUO2`
      - Class to set UO2 density to a constant input value
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/density/constantDensityMolybdenum`
      - Class to set Molybdenum density to a constant input value
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/density/densityModel`
      - Mother class for density models
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/density/densityConstant`
      - Model for constant density
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/density/densitySchumann1515Ti`
      - Class modelling the density evolution of 15-15Ti cladding material through 
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/density/constantDensityUPuO2`
      - Class modelling the constant density of (U,Pu)O2 MOX fuel
    * - :doc:`generated/offbeat/src/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/density/densityIAEAZy`
      - Class modelling Zy density from IAEA


Gap Gas Models
--------------

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/offbeat/src/offbeatLib/gapGasModel/gapTRISO
   generated/offbeat/src/offbeatLib/gapGasModel/gapGasTimeTabulated
   generated/offbeat/src/offbeatLib/gapGasModel/gapGasModel
   generated/offbeat/src/offbeatLib/gapGasModel/gapFRAPCON

.. list-table::
    :widths: 50 50

    * - :doc:`generated/offbeat/src/offbeatLib/gapGasModel/gapTRISO`
      - Model for gap gas composition, volume, temperature and pressure for TRISO
    * - :doc:`generated/offbeat/src/offbeatLib/gapGasModel/gapGasTimeTabulated`
      - Model for time-dependent gap gas composition, and pressure
    * - :doc:`generated/offbeat/src/offbeatLib/gapGasModel/gapGasModel`
      - As the name implies, the `none` gapGas class in OFFBEAT neglects the presence
    * - :doc:`generated/offbeat/src/offbeatLib/gapGasModel/gapFRAPCON`
      - This models interacts with the boundaries facing the gap to keep track of gap 


----------------
Function Objects
----------------

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/functionObjects/fieldDiffExtents/fieldDiffExtents
   generated/functionObjects/patchScalarFieldValue/patchScalarFieldValue
   generated/functionObjects/TBulk/TBulk
   generated/functionObjects/pressureDrop/pressureDrop
   generated/functionObjects/stopIfMaxFieldDiff/stopIfMaxFieldDiff
   generated/functionObjects/massFlow/massFlow
   generated/offbeat/src/functionObjects/centralVoidRadius/centralVoidRadius
   generated/offbeat/src/functionObjects/solidusRadius/solidusRadius
   generated/offbeat/src/functionObjects/gapGasComposition/gapGasComposition
   generated/offbeat/src/functionObjects/axialExpansion/axialExpansion
   generated/offbeat/src/functionObjects/fgp/fgp
   generated/offbeat/src/functionObjects/fgr/fgr
   generated/offbeat/src/functionObjects/rodPressure/rodPressure
   generated/offbeat/src/functionObjects/radialExpansion/radialExpansion
   generated/offbeat/src/functionObjects/radialProfileScalarField/radialProfileScalarField
   generated/offbeat/src/functionObjects/avgFieldBySlice/avgFieldBySlice
   generated/offbeat/src/functionObjects/columnarGrainRadius/columnarGrainRadius
   generated/offbeat/src/functionObjects/patchProfileScalarField/patchProfileScalarField
   generated/offbeat/src/functionObjects/patchRadius/patchRadius

.. list-table::
    :widths: 50 50

    * - :doc:`generated/functionObjects/fieldDiffExtents/fieldDiffExtents`
      - Calculates the spatial minimum and maximum extents of a difference between
    * - :doc:`generated/functionObjects/patchScalarFieldValue/patchScalarFieldValue`
      - Writes the values of an fvPatchScalarField
    * - :doc:`generated/functionObjects/TBulk/TBulk`
      - Calculates the total mass flux through a set of faceSets, faceZones or
    * - :doc:`generated/functionObjects/pressureDrop/pressureDrop`
      - Calculates the spatial minimum and maximum extents of a field
    * - :doc:`generated/functionObjects/stopIfMaxFieldDiff/stopIfMaxFieldDiff`
      - Stops the simulation when max(field1-field2) > 0 and prints out the
    * - :doc:`generated/functionObjects/massFlow/massFlow`
      - Calculates the total mass flux through a set of faceSets, faceZones or
    * - :doc:`generated/offbeat/src/functionObjects/centralVoidRadius/centralVoidRadius`
      - FunctionObject that prints central void radii vs time
    * - :doc:`generated/offbeat/src/functionObjects/solidusRadius/solidusRadius`
      - FunctionObject that prints solidus radius vs time in a specific cellZone
    * - :doc:`generated/offbeat/src/functionObjects/gapGasComposition/gapGasComposition`
      - FunctionObject that prints rod inner gas composition vs time
    * - :doc:`generated/offbeat/src/functionObjects/axialExpansion/axialExpansion`
      - FunctionObject that prints radius (assumed along x) of a patch accounting for
    * - :doc:`generated/offbeat/src/functionObjects/fgp/fgp`
      - FunctionObject that prints fission gas produced vs time
    * - :doc:`generated/offbeat/src/functionObjects/fgr/fgr`
      - FunctionObject that prints fgr% vs time
    * - :doc:`generated/offbeat/src/functionObjects/rodPressure/rodPressure`
      - FunctionObject that prints rod inner pressure vs time
    * - :doc:`generated/offbeat/src/functionObjects/radialExpansion/radialExpansion`
      - FunctionObject that prints radial expansion (assumed along x) of a patch 
    * - :doc:`generated/offbeat/src/functionObjects/radialProfileScalarField/radialProfileScalarField`
      - FunctionObject that prints rad profile of a scalar quantity at different axial 
    * - :doc:`generated/offbeat/src/functionObjects/avgFieldBySlice/avgFieldBySlice`
      - FunctionObject that prints avg of a scalar quantity at different axial 
    * - :doc:`generated/offbeat/src/functionObjects/columnarGrainRadius/columnarGrainRadius`
      - FunctionObject that prints columnar grain radius vs time
    * - :doc:`generated/offbeat/src/functionObjects/patchProfileScalarField/patchProfileScalarField`
      - FunctionObject that prints patch profile (for each face center) of a scalar 
    * - :doc:`generated/offbeat/src/functionObjects/patchRadius/patchRadius`
      - FunctionObject that prints radius (assumed along x) of a patch accounting for


-------------------
Boundary Conditions
-------------------



Neutronics
----------

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/fvPatchFields/neutronics/albedoSP3/albedoSP3FvPatchField

.. list-table::
    :widths: 50 50

    * - :doc:`generated/fvPatchFields/neutronics/albedoSP3/albedoSP3FvPatchField`
      - Albedo boundary condition for SP3 or diffusion calculations


Thermal-hydraulics
------------------

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/fvPatchFields/thermalHydraulics/blackBodyRadiation/blackBodyRadiationFvPatchScalarField
   generated/fvPatchFields/thermalHydraulics/mappedDarcyVelocity/mappedDarcyVelocityFvPatchVectorField
   generated/fvPatchFields/thermalHydraulics/NusseltThermalBaffle1D/NusseltThermalBaffle1DFvPatchScalarField
   generated/fvPatchFields/thermalHydraulics/velocityRundown/velocityRundownFvPatchVectorField
   generated/fvPatchFields/thermalHydraulics/timeFieldTable/timeFieldTableFvPatchScalarField
   generated/fvPatchFields/thermalHydraulics/fixedMassFlowRate/fixedMassFlowRateFvPatchVectorField
   generated/fvPatchFields/thermalHydraulics/parabolicVelocityPatchField/parabolicVelocityFvPatchVectorField

.. list-table::
    :widths: 50 50

    * - :doc:`generated/fvPatchFields/thermalHydraulics/blackBodyRadiation/blackBodyRadiationFvPatchScalarField`
      - This boundary condition provides a simplified black body radiation condition
    * - :doc:`generated/fvPatchFields/thermalHydraulics/mappedDarcyVelocity/mappedDarcyVelocityFvPatchVectorField`
      - Velocity inlet boundary condition either correcting the extrapolated
    * - :doc:`generated/fvPatchFields/thermalHydraulics/NusseltThermalBaffle1D/NusseltThermalBaffle1DFvPatchScalarField`
      - This BC models heat conduction through a thin 1-D baffle thermally coupling
    * - :doc:`generated/fvPatchFields/thermalHydraulics/velocityRundown/velocityRundownFvPatchVectorField`
      - Velocity inlet boundary condition either correcting the extrapolated
    * - :doc:`generated/fvPatchFields/thermalHydraulics/timeFieldTable/timeFieldTableFvPatchScalarField`
      - Boundary condition that allows to specify a scalarField over the boundary
    * - :doc:`generated/fvPatchFields/thermalHydraulics/fixedMassFlowRate/fixedMassFlowRateFvPatchVectorField`
      - Exactly the same as OpenFOAM's flowRateInletVelocity, but renormalizes the
    * - :doc:`generated/fvPatchFields/thermalHydraulics/parabolicVelocityPatchField/parabolicVelocityFvPatchVectorField`
      - Boundary condition specifies a parabolic velocity inlet profile


Thermo-mechanics
----------------

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/fvPatchFields/thermoMechanics/tractionDisplacement/gapContactFvPatchVectorField
   generated/fvPatchFields/thermoMechanics/tractionDisplacement/GFTractionDisplacementFvPatchVectorField

.. list-table::
    :widths: 50 50

    * - :doc:`generated/fvPatchFields/thermoMechanics/tractionDisplacement/gapContactFvPatchVectorField`
      - Fixed traction boundary condition for the standard linear elastic, fixed
    * - :doc:`generated/fvPatchFields/thermoMechanics/tractionDisplacement/GFTractionDisplacementFvPatchVectorField`
      - Fixed traction boundary condition for the standard linear elastic, fixed


Thermo-mechanics
----------------

.. toctree::
   :hidden:
   :maxdepth: 1

   generated/offbeat/src/offbeatLib/fvPatchFields/regionCoupledFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/zeroCurrent/hydrogenTransportFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/zeroCurrent/zeroCurrentActinidesRedistributionFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/zeroCurrent/porosityOutletFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/zeroCurrent/oxidePickupFractionFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/gaussianLaser/gaussianLaserFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/fixedDisplacement/fixedDisplacementFvPatchVectorField
   generated/offbeat/src/offbeatLib/fvPatchFields/solidDirectionMixed/solidDirectionMixedFvPatchVectorField
   generated/offbeat/src/offbeatLib/fvPatchFields/fixedTemperature/fixedTemperatureFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/fixedTemperature/timeDependentTemperatureFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/fixedTemperature/timeDependentAxialAzimuthalProfileFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/fixedTemperature/fmiAxialProfileFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/fixedTemperature/xzTemperatureProfileFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/fixedTemperature/timeDependentAzimuthalProfileFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/fixedTemperature/timeDependentAxialProfileFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/fixedTemperature/axialProfileFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/fixedTemperature/azimuthalProfileFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/radiativeConvectiveSink/radiativeConvectiveSinkFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/convectiveHTC/fmiAxialProfileHTCfvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/convectiveHTC/convectiveFromFluidFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/convectiveHTC/timeDependentAxialProfileHTCfvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/convectiveHTC/convectiveHTCFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/multiRegion/regionCoupledTemperatureFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/extAxialPatch/extAxialPatch
   generated/offbeat/src/offbeatLib/fvPatchFields/timeDependentAverageMapped/timeDependentAverageMappedFvPatchFieldsFwd
   generated/offbeat/src/offbeatLib/fvPatchFields/timeDependentAverageMapped/timeDependentAverageMappedFvPatchFields
   generated/offbeat/src/offbeatLib/fvPatchFields/timeDependentAverageMapped/timeDependentAverageMappedFvPatchField
   generated/offbeat/src/offbeatLib/fvPatchFields/extrapolatedFixedValue/extrapolatedFixedValueFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/tractionDisplacement/fixedTorqueFvPatchVectorField
   generated/offbeat/src/offbeatLib/fvPatchFields/tractionDisplacement/fixedPointLoadFvPatchVectorField
   generated/offbeat/src/offbeatLib/fvPatchFields/tractionDisplacement/topCladRingPressureFvPatchVectorField
   generated/offbeat/src/offbeatLib/fvPatchFields/tractionDisplacement/plenumSpringPressureFvPatchVectorField
   generated/offbeat/src/offbeatLib/fvPatchFields/tractionDisplacement/coolantPressureFvPatchVectorField
   generated/offbeat/src/offbeatLib/fvPatchFields/tractionDisplacement/contactFvPatchVectorField
   generated/offbeat/src/offbeatLib/fvPatchFields/tractionDisplacement/gapPressureFvPatchVectorField
   generated/offbeat/src/offbeatLib/fvPatchFields/tractionDisplacement/tractionDisplacementFvPatchVectorField
   generated/offbeat/src/offbeatLib/fvPatchFields/codedFixedTraction/codedTractionDisplacementFvPatchVectorField
   generated/offbeat/src/offbeatLib/fvPatchFields/codedFixedTraction/codedFixedGradientFvPatchFields
   generated/offbeat/src/offbeatLib/fvPatchFields/codedFixedTraction/codedFixedGradientFvPatchField
   generated/offbeat/src/offbeatLib/fvPatchFields/codedFixedTraction/codedFixedGradientFvPatchFieldsFwd
   generated/offbeat/src/offbeatLib/fvPatchFields/implicitContact/implicitGapContactFvPatchVectorField
   generated/offbeat/src/offbeatLib/fvPatchFields/temperatureCoupled/fuelRodGapFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/temperatureCoupled/temperatureCoupledFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/temperatureCoupled/resistiveGapFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/temperatureCoupled/trisoGapFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/fixedDisplacementZeroShear/fixedDisplacementZeroShearFvPatchVectorField
   generated/offbeat/src/offbeatLib/fvPatchFields/fixedDisplacementZeroShear/unilateralContactFvPatchVectorField
   generated/offbeat/src/offbeatLib/fvPatchFields/coolantChannel/NaCoolantChannelfvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/coolantChannel/coolantChannelfvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/fixedHeatFlux/fixedHeatFluxFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/fixedHeatFlux/axialProfileHeatFluxFvPatchScalarField
   generated/offbeat/src/offbeatLib/fvPatchFields/calculatedTotalField/calculatedTotalDisplacementFvPatchVectorField
   generated/offbeat/src/offbeatLib/fvPatchFields/plenumSpringBase/plenumSpringBase

.. list-table::
    :widths: 50 50

    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/regionCoupledFvPatchScalarField`
      - General region-coupled patchField
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/zeroCurrent/hydrogenTransportFvPatchScalarField`
      - Boundary condition to use with hydrogen transport
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/zeroCurrent/zeroCurrentActinidesRedistributionFvPatchScalarField`
      - Usage
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/zeroCurrent/porosityOutletFvPatchScalarField`
      - Boundary condition implemented for testing porosity outflow
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/zeroCurrent/oxidePickupFractionFvPatchScalarField`
      - Boundary condition to use with hydrogen transport
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/gaussianLaser/gaussianLaserFvPatchScalarField`
      - BOundary condition to prescribe a gaussian LASER intensity profile on a 
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/fixedDisplacement/fixedDisplacementFvPatchVectorField`
      - The `fixedDisplacement` boundary condition allows imposing a fixed displacement 
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/solidDirectionMixed/solidDirectionMixedFvPatchVectorField`
      - directionMixed with non-orthogonal correction for the diffusion term
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/fixedTemperature/fixedTemperatureFvPatchScalarField`
      - The `fixedTemperature` boundary condition allows to specify a constant fixed
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/fixedTemperature/timeDependentTemperatureFvPatchScalarField`
      - The `timeDependentTemperature` boundary condition allows to impose a time
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/fixedTemperature/timeDependentAxialAzimuthalProfileFvPatchScalarField`
      - Boundary condition to specify a time-dependent axial and azimuthal temperature profile
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/fixedTemperature/fmiAxialProfileFvPatchScalarField`
      - Boundary condition to impose a fixed axial temperature profile from an FMI
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/fixedTemperature/xzTemperatureProfileFvPatchScalarField`
      - The `xzTemperatureProfile` boundary condition allows to impose a temperature
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/fixedTemperature/timeDependentAzimuthalProfileFvPatchScalarField`
      - Boundary condition to impose a time-dependent azimuthal temperature profile
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/fixedTemperature/timeDependentAxialProfileFvPatchScalarField`
      - The `timeDependentAxialProfileT` boundary condition allows to impose a time
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/fixedTemperature/axialProfileFvPatchScalarField`
      - The `axialProfileT` boundary condition allows to impose a temperature axial
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/fixedTemperature/azimuthalProfileFvPatchScalarField`
      - Boundary condition to impose a temperature axial profile (constant in time)
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/radiativeConvectiveSink/radiativeConvectiveSinkFvPatchScalarField`
      - Heat sink boundary for solid temperature defined by radiative and
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/convectiveHTC/fmiAxialProfileHTCfvPatchScalarField`
      - Heat sink boundary defined by an axial profile of convective heat transfer
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/convectiveHTC/convectiveFromFluidFvPatchScalarField`
      - This class implements a mixed boundary condition which is meant
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/convectiveHTC/timeDependentAxialProfileHTCfvPatchScalarField`
      - The `timeDependentAxialProfileHTC` boundary condition allows to simulate an heat
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/convectiveHTC/convectiveHTCFvPatchScalarField`
      - The `convectiveHTC` boundary condition allows to simulate an heat sink boundary
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/multiRegion/regionCoupledTemperatureFvPatchScalarField`
      - turbulentTemperatureCoupledBaffleMixedFvPatchScalarField for OFFBEAT
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/extAxialPatch/extAxialPatch`
      - FMI function object to extract sampled list of values along a patch on the
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/timeDependentAverageMapped/timeDependentAverageMappedFvPatchFieldsFwd`
      - No description available
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/timeDependentAverageMapped/timeDependentAverageMappedFvPatchFields`
      - No description available
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/timeDependentAverageMapped/timeDependentAverageMappedFvPatchField`
      - mappedFixedValue boundary conditions where the average is defined using
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/extrapolatedFixedValue/extrapolatedFixedValueFvPatchScalarField`
      - Extrapolate internal field to the boundary as a fixedValue
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/tractionDisplacement/fixedTorqueFvPatchVectorField`
      - The `fixedTorque` boundary condition applies a fixed torque (in Nm) uniformly 
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/tractionDisplacement/fixedPointLoadFvPatchVectorField`
      - The `fixedPointLoad` boundary condition applies a fixed point load to a patch
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/tractionDisplacement/topCladRingPressureFvPatchVectorField`
      - The `topCladRingPressure` boundary condition is used for the top cladding 
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/tractionDisplacement/plenumSpringPressureFvPatchVectorField`
      - The `plenumSpringPressure` boundary condition replicates the counteracting 
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/tractionDisplacement/coolantPressureFvPatchVectorField`
      - The `coolantPressure` boundary condition is used to apply the coolant pressure 
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/tractionDisplacement/contactFvPatchVectorField`
      - The `contact` boundary condition uses the penalty method to handle
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/tractionDisplacement/gapPressureFvPatchVectorField`
      - The `gapPressure` boundary condition applies the gap pressure to the selected 
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/tractionDisplacement/tractionDisplacementFvPatchVectorField`
      - The `tractionDisplacement` boundary condition is used to apply a fixed or time-
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/codedFixedTraction/codedTractionDisplacementFvPatchVectorField`
      - Fixed traction boundary condition for the standard linear elastic, fixed
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/codedFixedTraction/codedFixedGradientFvPatchFields`
      - No description available
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/codedFixedTraction/codedFixedGradientFvPatchField`
      - Constructs on-the-fly a new boundary condition (derived from
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/codedFixedTraction/codedFixedGradientFvPatchFieldsFwd`
      - No description available
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/implicitContact/implicitGapContactFvPatchVectorField`
      - The `implicitGapContact` boundary condition is an extension of the standard 
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/temperatureCoupled/fuelRodGapFvPatchScalarField`
      - Coupled boundary condition for modeling the presence of an evolving gap
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/temperatureCoupled/temperatureCoupledFvPatchScalarField`
      - Coupled boundary condition that enforces continuity for the temperature field
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/temperatureCoupled/resistiveGapFvPatchScalarField`
      - Coupled boundary condition for modeling the presence of a fixed gap (or contact)
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/temperatureCoupled/trisoGapFvPatchScalarField`
      - Region coupled patchField for temperature fields in solid heat conduction in
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/fixedDisplacementZeroShear/fixedDisplacementZeroShearFvPatchVectorField`
      - The `fixedDisplacementZeroShear` boundary condition fixes the displacement 
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/fixedDisplacementZeroShear/unilateralContactFvPatchVectorField`
      - The `unilateralContact` boundary condition enforces a contact-like behavior 
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/coolantChannel/NaCoolantChannelfvPatchScalarField`
      - Coolant subchannel module to provide thermal boundary condition for Na 
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/coolantChannel/coolantChannelfvPatchScalarField`
      - Coolant subchannel module to provide thermal boundary condition
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/fixedHeatFlux/fixedHeatFluxFvPatchScalarField`
      - Fixed heat flux boundary condition
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/fixedHeatFlux/axialProfileHeatFluxFvPatchScalarField`
      - Fixed heat flux boundary condition
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/calculatedTotalField/calculatedTotalDisplacementFvPatchVectorField`
      - This BC is created to be the default BC when the displacement field D is 
    * - :doc:`generated/offbeat/src/offbeatLib/fvPatchFields/plenumSpringBase/plenumSpringBase`
      - Base class with common functinality for plenumSpring patch fields


