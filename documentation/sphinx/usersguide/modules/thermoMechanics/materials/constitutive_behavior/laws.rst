Constitutive laws
=================

This section lists the constitutive laws available in OFFBEAT and in the
thermo-mechanics module used within multi-physics GeN-Foam simulations.

Constitutive laws define the stress–strain response of a material,
including elasticity, plasticity, creep, and their possible
combinations.

The constitutive law is selected inside the material definition in the
``materials`` section of ``solverDict``.

The model is selected using the ``type`` keyword inside the
``rheology`` subdictionary.

.. code-block:: cpp

   materials
   {
       cladding
       {
           type Zircaloy;

           rheology
           {
               type misesPlasticCreep;
           }
       }
   }

Depending on the selected law, additional submodels may be required,
such as yield stress models or creep correlations.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following
constitutive laws:

.. toctree::
   :maxdepth: 1
   :caption: Constitutive laws
   :glob:

   ../../../../../../cppapi/generated/offbeatLib/rheology/constitutiveLaws/*