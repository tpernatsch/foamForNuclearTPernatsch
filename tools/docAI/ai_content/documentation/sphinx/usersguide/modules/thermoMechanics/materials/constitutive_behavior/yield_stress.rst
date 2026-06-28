Yield stress models
===================

Yield stress models provide the yield strength used by plasticity-based
constitutive laws.

They typically define how the yield stress depends on temperature,
irradiation, accumulated plastic strain, strain rate, or other state
variables, depending on the selected model.

In OFFBEAT and in the thermo-mechanics module used within multi-physics
GeN-Foam simulations, yield stress models are defined inside the
constitutive law configuration in the ``materials`` section of
``solverDict``.

The model is selected using the ``type`` keyword inside the
``yieldStress`` subdictionary of the chosen constitutive law.

.. code-block:: cpp

   materials
   {
       cladding
       {
           type Zircaloy;

           rheology
           {
               type misesPlasticity;

               yieldStress
               {
                   type Fraptran;
               }
           }
       }
   }

Depending on the selected model, the yield stress may evolve with
temperature, irradiation conditions, or accumulated plastic strain. It
also directly affects the onset of plastic deformation in the material.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following yield
stress models:

.. toctree::
   :maxdepth: 1
   :caption: Yield stress models
   :glob:

   ../../../../../../cppapi/generated/offbeatLib/rheology/constitutiveLaws/yieldStressModels/*