Young modulus
=============

The Young modulus defines the elastic stiffness of a material and is a key
mechanical property used by the mechanical solver in OFFBEAT and in the
thermo-mechanics module used within multi-physics GeN-Foam simulations.

The Young modulus model is defined inside the material subdictionary
of the ``materials`` section in ``solverDict``.

The model is selected using the ``type`` keyword inside the
``YoungModulus`` property subdictionary.

.. code-block:: cpp

   materials
   {
       fuel
       {
           type UO2;

           YoungModulus
           {
               type UO2Matpro;
           }
       }
   }

The selected model provides the temperature-dependent elastic modulus
used by the mechanical solver.

-------------------------------------------------------------------------------

.. rubric:: Constant property shortcut

For materials with constant properties, the Young modulus may be
specified directly without defining a property subdictionary.

.. code-block:: cpp

   materials
   {
       cladding
       {
           type constant;

           E [1 -1 -2 0 0 0 0] 1.0e11;
       }
   }

In this case, the Young modulus is treated as constant and no additional
model needs to be selected.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following Young
modulus models:

.. toctree::
   :maxdepth: 1
   :caption: Young modulus models
   :glob:

   ../../../../../../cppapi/generated/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/YoungModulus/*