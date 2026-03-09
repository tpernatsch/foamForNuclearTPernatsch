Creep models
============

Creep models provide the creep strain rate (or equivalent creep
contribution) used by constitutive laws that include time-dependent
inelastic deformation.

Creep represents the gradual accumulation of deformation under sustained
stress, typically driven by temperature and irradiation. In nuclear fuel
and structural materials, creep can play a major role in long-term
deformation and stress redistribution.

In OFFBEAT and in the thermo-mechanics module used within multi-physics
GeN-Foam simulations, creep models are defined inside the constitutive
law configuration in the ``materials`` section of ``solverDict``.

The model is selected using the ``type`` keyword inside the
``creep`` subdictionary of the chosen constitutive law.

.. code-block:: cpp

   materials
   {
       cladding
       {
           type Zircaloy;

           rheology
           {
               type misesPlasticCreep;

               creep
               {
                   type ZircaloyLimback;
               }
           }
       }
   }

Depending on the constitutive law, the creep strain rate may depend on
temperature, stress, irradiation conditions, or other material state
variables.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following creep
models:

.. toctree::
   :maxdepth: 1
   :caption: Creep models
   :glob:

   ../../../../../../cppapi/generated/offbeatLib/rheology/constitutiveLaws/creepModels/*