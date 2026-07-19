Heat source
===========

The ``heatSource`` model handles the definition of the volumetric heat
generation in the thermo-mechanics module (OFFBEAT / GeN-Foam solid region).

The power density field is named ``Q`` and is expressed in ``W/m^3``.

-------------------------------------------------------------------------------

.. rubric:: Selection in ``solverDict``

The heat source model is selected in ``solverDict`` via the ``heatSource``
dictionary (located in ``constant/`` for OFFBEAT, or in
``constant/<regionName>/`` for GeN-Foam):

.. code-block:: cpp

   heatSource
   {
       type constantLhgr;
       // additional options (if any)
   }

If the ``heatSource`` dictionary is not present in ``solverDict``, the heat
source is neglected and the field ``Q`` is not created.

-------------------------------------------------------------------------------

.. note::

   Traditional 1D fuel performance codes typically require the radially
   averaged linear heat generation rate (``lhgr``) as input, often together
   with an axial profile.

   For 3D simulations with arbitrary geometries and unstructured meshes,
   defining the heat source field is less straightforward.

   When the power density field is not symmetric, the most common approach
   is to couple the thermo-mechanics module with a neutronics or multiphysics
   solver that provides the full 3D power distribution. In that case, the
   heat source field can be read directly from a previous simulation.

   Alternatively, OpenFOAM utilities such as ``topoSet`` can be used to
   construct spatially varying heat source fields.

   For simulations where the heat source is assumed uniform along the
   azimuthal direction, the ``lhgr``-based models listed below can be used.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following heat source
models:

.. toctree::
   :maxdepth: 1
   :caption: Heat source models
   :glob:

   ../../../../../cppapi/generated/offbeatLib/heatSource/*