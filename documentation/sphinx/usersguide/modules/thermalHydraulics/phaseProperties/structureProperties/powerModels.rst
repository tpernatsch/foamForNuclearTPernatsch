.. _modules_thermalHydraulics_porousMdeium_powerModels:

Power models
^^^^^^^^^^^^

Power models are  models used to simulate the behavior  of sub-scale structures that produce power. Main entries are repoted below.

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Entry
     - Mandatory
     - Description
   * - ``type``
     - Y
     - Run-type selectable model (see below for a list of available models)
   * - ``volumeFraction``
     - N
     - Used to compute the volumetric heat capacity of the structure as volumeFraction*rho*Cp.  Defaults to the overall structure volumeFraction the start of the cellZone subDict
   * - ``volumetricArea``
     - N/Y
     - Structure surface area per unit volume. Normally needed. Not needed e.g., for nuclearFuelPin and heatedPin power models  since it can be derived by the volumeFraction
   * - ``powerOffCriterionModel``
     - N
     - See :ref:`powerOffCriterionModel <modules_thermalHydraulics_porousMedium_powerOffCriterionModels>`

Other entries depends on the chosen power model. 

An example of powerModel sub-dictionary is reported below:

.. code :: cpp

    powerModel
    {
        type            heatedPin;
        innerRadius     0;
        outerRadius     0.003;
        k               20;
        meshSize        6;
        T               653.15;
        Cp              500;
        rho             7700;

        powerOffCriterionModel
        {
            type            timeThreshold;
            time            12.475;
        }
    }


Available models include:


