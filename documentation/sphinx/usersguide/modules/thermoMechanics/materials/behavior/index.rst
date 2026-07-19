Behavioral models
=================

Behavioral models describe **material evolution driven by operating
conditions** such as temperature, stress, irradiation, or burnup.

Some behavioral models contribute an additional strain component
(e.g. swelling or densification), while others act through state
variables or failure criteria and do not directly modify the strain
field.

In OFFBEAT and in the thermo-mechanics module used within multi-physics
GeN-Foam simulations, behavioral models are defined **per material**
inside the ``materials`` dictionary of ``solverDict``.

Typical examples include fuel swelling, densification, relocation,
phase transitions, and failure criteria.

.. note::

   Not all materials support all behavioral models.
   The set of available behavioral models depends on the selected
   material type, and some materials may not use any behavioral
   model at all.

-------------------------------------------------------------------------------

The following pages document the available behavioral models.

.. toctree::
   :maxdepth: 1
   :caption: Behavioral models

   densification
   relocation
   swelling
   failure
   phase_transition