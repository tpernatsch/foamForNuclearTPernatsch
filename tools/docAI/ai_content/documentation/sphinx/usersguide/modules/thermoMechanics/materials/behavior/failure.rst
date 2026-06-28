Failure models
==============

Failure models define criteria for material failure, loss of integrity,
or detection of unphysical conditions during a simulation.

These models typically evaluate stress, strain, temperature, or other
state variables to determine whether a material has reached a failure
limit. Depending on the model, failure may correspond to phenomena
such as melting, plastic instability, excessive strain, or other
integrity criteria.

In OFFBEAT and in the thermo-mechanics module used within multi-physics
GeN-Foam simulations, failure models are defined inside the material
subdictionary of the ``materials`` section in ``solverDict``.

The model is selected using the ``type`` keyword inside the
``failure`` subdictionary.

.. code-block:: cpp

   materials
   {
       cladding
       {
           type Zircaloy;

           failure
           {
               type ZircaloyOverstrainBison;
           }
       }
   }

Depending on the selected model, failure conditions may trigger
warnings, terminate the simulation, or provide diagnostic information
indicating that a material limit has been exceeded.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following failure
models:

.. toctree::
   :maxdepth: 1
   :caption: Failure models
   :glob:

   ../../../../../../cppapi/generated/offbeatLib/materials/materialModel/behavioralModels/failureModels/*