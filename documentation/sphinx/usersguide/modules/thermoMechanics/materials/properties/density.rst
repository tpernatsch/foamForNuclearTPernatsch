Density
=======

Density defines the mass per unit volume of a material and is a
fundamental thermo-mechanical property used in OFFBEAT and in the
thermo-mechanics module used within multi-physics GeN-Foam simulations.

The density model is defined inside the material subdictionary
of the ``materials`` section in ``solverDict``.

The model is selected using the ``type`` keyword inside the
``density`` property subdictionary.

.. code-block:: cpp

   materials
   {
       fuel
       {
           type UO2;

           density
           {
               type UO2Constant;
           }
       }
   }

Depending on the selected model, density may be treated as constant or
computed using empirical correlations.

-------------------------------------------------------------------------------

.. rubric:: Constant property shortcut

For materials with constant properties, the density may be specified
directly without defining a property subdictionary.

.. code-block:: cpp

   materials
   {
       cladding
       {
           type constant;

           rho rho [1 -3 0 0 0 0 0] 6500;
       }
   }

In this case, the density is treated as constant and no additional
model needs to be selected.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following density
models:

.. toctree::
   :maxdepth: 1
   :caption: Density models
   :glob:

   ../../../../../../cppapi/generated/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/density/*