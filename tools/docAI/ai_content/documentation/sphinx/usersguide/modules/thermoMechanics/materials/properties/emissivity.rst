Emissivity
==========

Emissivity defines the efficiency of a material surface in emitting
thermal radiation. It is used by radiative heat transfer models in
OFFBEAT and in the thermo-mechanics module used within multi-physics
GeN-Foam simulations.

The emissivity model is defined inside the ``materials`` subdictionary
of ``solverDict``.

The model is selected using the ``type`` keyword inside the
``emissivity`` property subdictionary.

.. code-block:: cpp

   materials
   {
       fuel
       {
           type UO2;

           emissivity
           {
               type UO2Relap;
           }
       }
   }

Depending on the selected model, emissivity may be treated as constant
or computed using empirical correlations, depending on temperature and
material type.

-------------------------------------------------------------------------------

.. rubric:: Constant property shortcut

For materials with constant properties, emissivity may be specified
directly, without defining a property subdictionary.

.. code-block:: cpp

   materials
   {
       cladding
       {
           type constant;

           emissivity emissivity [0 0 0 0 0 0 0] 0.8;
       }
   }

In this case, emissivity is treated as constant, and no additional
model needs to be selected.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following emissivity
models:

.. toctree::
   :maxdepth: 1
   :caption: Emissivity models
   :glob:

   ../../../../../../cppapi/generated/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/emissivity/*