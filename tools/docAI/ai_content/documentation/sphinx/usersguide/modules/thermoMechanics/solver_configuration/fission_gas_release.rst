Fission gas release
===================

The ``fissionGasRelease`` model handles the generation, diffusion, and
release of gaseous fission products in the thermo-mechanics module
(OFFBEAT / GeN-Foam solid region).

-------------------------------------------------------------------------------

.. rubric:: Selection in ``solverDict``

The fission gas release model is selected in ``solverDict`` via the ``fgr``
dictionary (located in ``constant/`` for OFFBEAT, or in
``constant/<regionName>/`` for GeN-Foam):

.. code-block:: cpp

   fgr
   {
       type SCIANTIX;
       // additional options (if any)
   }

If the ``fgr`` dictionary is not present in ``solverDict``, fission gas
release is neglected, and the gap gas composition remains as defined in
the ``0/gapGas`` file.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following fission gas
release models:

.. toctree::
   :maxdepth: 1
   :caption: Fission gas release models
   :glob:

   ../../../../../cppapi/generated/offbeatLib/fissionGasRelease/*