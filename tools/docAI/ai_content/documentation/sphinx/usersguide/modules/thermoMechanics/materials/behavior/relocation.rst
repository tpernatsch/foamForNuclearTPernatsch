Relocation
==========

Relocation models describe the outward displacement of fuel material,
typically associated with cracking and fragmentation during the early
stages of irradiation.

This phenomenon results from the formation of radial cracks in the fuel.
These cracks allow fragments of the pellet to move outward and partially
close the fuel–cladding gap.

In OFFBEAT and in the thermo-mechanics module used within multi-physics
GeN-Foam simulations, relocation models are defined inside the ``materials``
subdictionary of the ``solverDict``.

The model is selected using the ``type`` keyword inside the
``relocation`` subdictionary.

.. code-block:: cpp

   materials
   {
       fuel
       {
           type UO2;

           relocation
           {
               type UO2Frapcon;
           }
       }
   }

The resulting relocation strain contributes to the overall deformation
of the fuel pellet and affects the evolution of the fuel–cladding gap.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following relocation
models:

.. toctree::
   :maxdepth: 1
   :caption: Relocation models
   :glob:

   ../../../../../../cppapi/generated/offbeatLib/materials/materialModel/behavioralModels/relocation/*