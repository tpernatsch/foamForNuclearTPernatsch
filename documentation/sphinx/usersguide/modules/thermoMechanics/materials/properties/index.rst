Thermo-mechanical properties
============================

This section documents the **thermo-mechanical material properties**
available in OFFBEAT and in the thermo-mechanics module used within
multi-physics GeN-Foam simulations.

These properties are used by the thermal and mechanical solvers and are
defined **per material** inside the ``materials`` dictionary of
``solverDict``.

Thermo-mechanical properties represent intrinsic material quantities or
correlations such as:

* thermal conductivity,
* heat capacity,
* density,
* emissivity,
* thermal expansion,
* elastic constants.

Each property is configured through a dedicated property subdictionary
inside the material definition. The model used to evaluate the property
is selected via the ``type`` keyword.

.. code-block:: cpp

   materials
   {
       fuel
       {
           type UO2;

           conductivity
           {
               type UO2Nfir;
           }

           Cp
           {
               type UO2Cp;
           }
       }
   }

The selected model may depend on temperature, burnup, porosity,
composition, or other state variables depending on the specific
correlation.

-------------------------------------------------------------------------------

The following pages describe the available thermo-mechanical
property models.

.. toctree::
   :maxdepth: 1
   :caption: Contents

   thermal_conductivity
   heat_capacity
   density
   emissivity
   thermal_expansion
   young_modulus
   poisson_ratio