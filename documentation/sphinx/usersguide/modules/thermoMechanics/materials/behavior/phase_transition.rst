Phase transition
================

Phase transition models describe changes in material phase induced by
temperature or transient operating conditions.

Such models are typically used for materials that may undergo structural
transformations during operation. For example, zirconium alloys may
experience phase changes at high temperatures, which can influence the
material properties and mechanical response.

In OFFBEAT and in the thermo-mechanics module used within multi-physics
GeN-Foam simulations, phase transition models are defined inside the
material subdictionary of the ``materials`` section in ``solverDict``.

The model is selected using the ``type`` keyword inside the
``phaseTransition`` subdictionary.

.. code-block:: cpp

   materials
   {
       cladding
       {
           type Zircaloy;

           phaseTransition
           {
               type ZircaloyDynamic;
           }
       }
   }

Depending on the selected model, the phase state of the material may
affect the evaluation of other material properties or mechanical
behavior.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following phase
transition models:

.. toctree::
   :maxdepth: 1
   :caption: Phase transition models
   :glob:

   ../../../../../../cppapi/generated/offbeatLib/materials/materialModel/behavioralModels/phaseTransition/*