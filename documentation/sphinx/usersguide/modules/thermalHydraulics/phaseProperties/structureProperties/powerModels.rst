.. _modules_thermalHydraulics_porousMdeium_powerModels:

Power models
^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

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

Other entries depends on the chosen power model. Available models include:
