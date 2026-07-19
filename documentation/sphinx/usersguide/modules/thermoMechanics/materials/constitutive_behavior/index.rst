Constitutive mechanical behavior
=================

Constitutive laws define the mechanical response of a material, i.e.
how stress is computed from the strain history and internal variables
(elasticity, plasticity, creep, and their combinations).

In OFFBEAT and in the thermo-mechanics module used within multi-physics
GeN-Foam simulations, the constitutive law is defined **per material**
inside the ``materials`` dictionary of ``solverDict``.

Options that are related to the general rheological framework
(e.g. plane stress/strain assumptions and optional thermal stress
settings) are defined in the ``rheology`` dictionary of
``solverDict``. These options are documented in
:doc:`../../solver_configuration/rheology`.

-------------------------------------------------------------------------------

The following pages describe the available constitutive models and
their associated submodels.

.. toctree::
   :maxdepth: 1
   :caption: Constitutive models

   laws
   yield_stress
   creep