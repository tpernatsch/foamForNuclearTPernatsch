Densification
=============

Densification models describe the reduction of material volume due to
pore closure, typically occurring during the early stages of irradiation.

This phenomenon is particularly relevant for oxide fuels, where the
as-fabricated porosity decreases as irradiation progresses, resulting
in a reduction of the fuel volume and a corresponding strain contribution.

In OFFBEAT and in the thermo-mechanics module used within multi-physics
GeN-Foam simulations, densification models are defined inside the
material subdictionary of the ``materials`` section in ``solverDict``.

The model is selected using the ``type`` keyword inside the
``densification`` subdictionary.

.. code-block:: cpp

   materials
   {
       fuel
       {
           type UO2;

           densification
           {
               type UO2Frapcon;
           }
       }
   }

The resulting densification strain contributes to the overall strain
state of the material.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following
densification models:

.. toctree::
   :maxdepth: 1
   :caption: Densification models
   :glob:

   ../../../../../../cppapi/generated/offbeatLib/materials/materialModel/behavioralModels/densification/*