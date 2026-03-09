Swelling
========

Swelling models describe volumetric changes induced by irradiation,
temperature, or microstructural evolution.

In nuclear fuels and structural materials, swelling is typically driven
by fission-product accumulation, gas bubble formation, or irradiation-
induced defect evolution. Depending on the material and model, swelling
may contribute an additional strain component to the mechanical response.

In OFFBEAT and in the thermo-mechanics module used within multi-physics
GeN-Foam simulations, swelling models are defined inside the material
subdictionary of the ``materials`` section in ``solverDict``.

The model is selected using the ``type`` keyword inside the
``swelling`` subdictionary.

.. code-block:: cpp

   materials
   {
       fuel
       {
           type UO2;

           swelling
           {
               type UO2Frapcon;
           }
       }
   }

The resulting swelling strain contributes to the
overall deformation of the material and can significantly affect gap
closure, stresses, and structural integrity.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following swelling
models:

.. toctree::
   :maxdepth: 1
   :caption: Swelling models
   :glob:

   ../../../../../../cppapi/generated/offbeatLib/materials/materialModel/behavioralModels/swelling/*