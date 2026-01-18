.. _modules_thermalHydraulics_porousMedium_twoPhaseSpecific:

Sub-dictionaries specific to two-phase simulations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Two-phase flow simulations require the use of few additional keywords and dictionaries compared to single-phase simulations.
- The keywords ``fluid1`` and ``fluid2`` are used to define arbitrary names for fluids 1 and 2. Fields in the time folders have  these in their names (after the dot) and everything else must be named accordingly in the phaseProperties dictionary.
- Two sub-dictionaries, one for each fluid. These dictionaries should always be named fluid1NameProperties and fluid2NameProperties, with the names fluid1 and fluid2 defined above.

Each of the two subdictionaries should contain the following keywords

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Entry
     - Mandatory
     - Description
   * - stateOfMatter
     - Y/N
     - Supported entries are either "gas" or "liquid". This entry is used by some drag, heat or mass              transfer models as knowledge of which phase is gaseous(/vapourous) and which phase is liquid is necessary sometimes. This entry is not mandatory per se, but might be required by specific choices of models.
   * - ``residualAlpha``
     - N
     - Residual fluid volume fraction used to stabilize quations when/if its volumeFraction tends to 0. Defaults to 1e-9
   * - ``thermoResidualAlpha``
     - N
     - Residual fluid volume fraction  below which the fluid temperature is not obtained from the enthalpy equation solution, but is set to  the fluid1-fluid2 interfacial temperature. This is applied on a cell-by-cell basis (e.g. only in those cells whose volumeFraction of this fluid is below thermoResidualAlpha). This is meant to be used to stabilize the temperature field of a fluid being produced during phase change (either boiling or condensation). In fact, small inaccuracies in the calculation of the absolute enthalpy of a phase with a small volume fraction will translate into large inaccuracies in the  temperature. Note that thermoResidualAlpha is expressed relatively to the available volume for fluid flow, thus accounting for possible structures. Defaults to 0.
   * - ``writeRestartFields``
     - N
     - If true, write additional fields to disk that quickly allow residuals to converge to their pre-restart values if restarting the simulation from a certain time-step. Defaults to false
   * - ``dispersedDiameterModel``
     - Y
     - Model used to calculate the diameter of the phase bubbles or droplets when the phase is dispersed. Options available are the templated ``constant``, ``byRegime``, or the specific models available in foamForNuclear (see :ref:`fluid diameter models <modules_thermalHydraulics_porousMedium_fluidDiameterModels>`)

An example of dictionary is reported below:

.. code :: cpp

    fluid1 "liquid";
    fluid2 "vapour";

    liquidProperties
    {
        stateOfMatter   liquid;
        residualAlpha 1e-9;
        dispersedDiameterModel
        {
            type    constant;
            value   0.01;
        }
    }

    vapourProperties
    {
        stateOfMatter       gas;
        thermoResidualAlpha 0.1;
        dispersedDiameterModel
        {
            type    constant;
            value   0.01;
        }
    }


