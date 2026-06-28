Poisson ratio
=============

The Poisson ratio defines the ratio between transverse and axial strain
under uniaxial loading. It is a basic elastic material property used in
OFFBEAT and in the thermo-mechanics module within multi-physics
GeN-Foam simulations.

The Poisson ratio model is defined inside the ``materials`` subdictionary
of ``solverDict``.

The model is selected using the ``type`` keyword inside the
``PoissonRatio`` property subdictionary.

.. code-block:: cpp

   materials
   {
       fuel
       {
           type UO2;

           PoissonRatio
           {
               type UO2Constant;
           }
       }
   }

Depending on the selected model, the Poisson ratio may be treated as
constant or defined through empirical correlations.

-------------------------------------------------------------------------------

.. rubric:: Constant property shortcut

For materials with constant properties, the Poisson ratio may be
specified directly without defining a property subdictionary.

.. code-block:: cpp

   materials
   {
       cladding
       {
           type constant;

           nu nu [0 0 0 0 0 0 0] 0.3;
       }
   }

In this case, the Poisson ratio is treated as constant, and no additional
model needs to be selected.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following Poisson
ratio models:

.. toctree::
   :maxdepth: 1
   :caption: Poisson ratio models
   :glob:

   ../../../../../../cppapi/generated/offbeatLib/materials/materialModel/thermoMechanicalPropertiesModels/PoissonRatio/*