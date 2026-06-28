Axial slice mapper
==================

The ``sliceMapper`` model defines axial slices that allow 1-D,
slice-based operations to be performed on a 3-D mesh in the
thermo-mechanics module (OFFBEAT / GeN-Foam solid region).

Although OFFBEAT supports arbitrary geometries and unstructured meshes,
some behavioural models in fuel performance originate from traditional
1.5-D approaches. These models may require quantities averaged over
axial slices, such as burnup or heat source density.

For example:

* certain relocation correlations
* the ``Lassmann`` burnup model
* axial force balance when using the ``modifiedPlaneStrain`` option in
  ``rheologyOptions``

In these situations, slice-averaged quantities must be computed in a
way that remains compatible with traditional fuel performance approaches.

However, OpenFOAM does not natively provide a concept of axial slices,
and such a concept cannot be directly generalized to arbitrary 3-D
domains. The ``sliceMapper`` class addresses this by defining a mapping
between the 3-D mesh cells and a set of virtual 1-D axial slices.

-------------------------------------------------------------------------------

.. rubric:: Selection in ``solverDict``

The axial slice mapper is selected in ``solverDict`` via the
``sliceMapper`` dictionary (located in ``constant/`` for OFFBEAT, or in
``constant/<regionName>/`` for GeN-Foam):

.. code-block:: cpp

   sliceMapper
   {
       type byMaterial;
       // additional options (if any)
   }

If the ``sliceMapper`` dictionary is not present in ``solverDict``, no
slice mapping is created. Some models that rely on slice-averaged
quantities may therefore not function correctly.

-------------------------------------------------------------------------------

.. warning::

   Axial slices are **not physical subdivisions of the mesh**.

   The ``sliceMapper`` instead constructs a *virtual mapping* by storing
   the indices of the cells belonging to each slice. Other models can
   then perform slice-based operations using this mapping as if the
   slices were explicitly present in the geometry.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following slice mapper
models:

.. toctree::
   :maxdepth: 1
   :caption: Slice mapper models
   :glob:

   ../../../../../cppapi/generated/offbeatLib/sliceMapper/*