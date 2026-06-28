Heat capacity
=============

Heat capacity defines the amount of energy required to change the
temperature of a material. It is a fundamental thermo-mechanical
property used in OFFBEAT and in the thermo-mechanics module within
multi-physics GeN-Foam simulations.

The heat capacity model is defined inside the material subdictionary
of the ``materials`` section in ``solverDict``.

The model is selected using the ``type`` keyword inside the
``heatCapacity`` property subdictionary.

.. code-block:: cpp

   materials
   {
       fuel
       {
           type UO2;

           heatCapacity
           {
               type UO2Matpro;
           }
       }
   }

The selected model provides the heat capacity as a function of
temperature. Depending on the specific correlation, it may also
depend on composition or other material state variables.

-------------------------------------------------------------------------------

.. rubric:: Constant property shortcut

For materials with constant properties, the heat capacity may be
specified directly without defining a property subdictionary.

.. code-block:: cpp

   materials
   {
       cladding
       {
           type constant;

           Cp Cp [0 2 -2 -1 0 0 0] 330;
       }
   }

In this case, the heat capacity is treated as constant, and no
additional model needs to be selected.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following heat
capacity models:

.. toctree::
   :maxdepth: 1
   :caption: Heat capacity models
   :glob:

   ../../../../../../cppapi/generated/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/heatCapacity/*