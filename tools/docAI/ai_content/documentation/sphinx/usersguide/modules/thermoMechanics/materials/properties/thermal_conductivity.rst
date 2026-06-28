Thermal conductivity
====================

Thermal conductivity defines a material's ability to conduct heat. It is a core material property used by the thermal solver in OFFBEAT and by the thermo-mechanics module within multi-physics GeN-Foam simulations.

The thermal conductivity model is defined inside the ``materials`` subdictionary of ``solverDict``.

The model is selected using the ``type`` keyword inside the ``conductivity`` property subdictionary.

.. code-block:: cpp

   materials
   {
       fuel
       {
           type UO2;

           conductivity
           {
               type UO2Nfir;
           }
       }
   }

The selected model provides the thermal conductivity as a function of temperature. Depending on the specific correlation, it may also depend on additional state variables, such as burnup, porosity, or composition.

-------------------------------------------------------------------------------

.. rubric:: Constant property shortcut

For materials with constant properties, the conductivity may be specified directly without defining a property subdictionary.

.. code-block:: cpp

   materials
   {
       cladding
       {
           type constant;

           k k [1 -1 -3 -1 0 0 0] 16;
       }
   }

In this case, the conductivity is treated as constant, and no additional model needs to be selected.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following thermal conductivity models:

.. toctree::
   :maxdepth: 1
   :caption: Thermal conductivity models
   :glob:

   ../../../../../../cppapi/generated/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/conductivity/*