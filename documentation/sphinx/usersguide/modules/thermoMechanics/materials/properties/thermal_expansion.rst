Thermal expansion
=================

Thermal expansion defines the change in material dimensions as a function
of temperature and is a key thermo-mechanical property used in OFFBEAT
and in the thermo-mechanics module used within multi-physics GeN-Foam
simulations.

The thermal expansion model is defined inside the material subdictionary
of the ``materials`` section in ``solverDict``.

The model is selected using the ``type`` keyword inside the
``thermalExpansion`` property subdictionary.

.. code-block:: cpp

   materials
   {
       fuel
       {
           type UO2;

           alphaT
           {
               type UO2Relap;
               Tref 300;
           }
       }
   }

The selected model provides the thermal strain contribution induced by
temperature changes.

-------------------------------------------------------------------------------

.. rubric:: Constant property shortcut

For materials with constant properties, the thermal expansion coefficient
may be specified directly without defining a property subdictionary.

.. code-block:: cpp

   materials
   {
       cladding
       {
           type constant;

           alpha alpha [0 0 0 -1 0 0 0] 1.0e-5;
           Tref Tref [0 0 0 1 0 0 0] 300;
       }
   }

In this case, the thermal expansion coefficient is treated as constant
and no additional model needs to be selected.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following thermal
expansion models:

.. toctree::
   :maxdepth: 1
   :caption: Thermal expansion models
   :glob:

   ../../../../../../cppapi/generated/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/thermalExpansion/*