.. _modules_thermalHydraulics_porousMedium_physicsModelsTwoPhase:

The *physicsModels* sub-dictionary for two-phase flows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The *physicsModels* sub-dictionary is used in two-phase simulations to define the models that describe the interaction between fluid1 and structure, fluid2 and structure, as well as fluid1 and fluid 2. Several types of models can be defined:



 drag models via the ``dragModels`` sub-dictionary and heat transfer models via the ``heatTransferModels`` sub-dictionary. Each of these two sub-dictionaries consists of a series of other sub-dictionaries whose properties are applied to the cell zones that have their name in the sub-dictionary keys. For both the drag models and the heat transfer models, two options exist:

- Use as type a single model chosen among those available (:ref:`drag models <modules_thermalHydraulics_porousMdeium_FSdrag>`, :ref:`heat transfer models <modules_thermalHydraulics_porousMedium_FSHeatTransferCoefficientModels>`)
- Use ``byRegime`` as  type, and then chose one of the available models (:ref:`drag models <modules_thermalHydraulics_porousMdeium_FSdrag>`, :ref:`heat transfer models <modules_thermalHydraulics_porousMedium_FSHeatTransferCoefficientModels>`) for every regime. This requires having set a :ref:`regime map <modules_thermalHydraulics_porousMedium_regimeMapModels>`.



An example of dictionary is reported below. 

.. code :: cpp
   
   physicsModels
   {
      dragModels
      {
         "diagrid:axial



The main run-time selctable models include:

.. toctree::
   :maxdepth: 1

   FFHeatTransferCoefficientModels
   FSHeatTransferCoefficientModels
   FFdrag
   FSdrag
   multipliersDrag
   phaseChangeModels
   contactPartitionModels
   dispersionModels
   fluidDiameterModels
   interfacialAreaModels
   virtualMassCoefficientModels
   latentHeatModels
   saturationModels



