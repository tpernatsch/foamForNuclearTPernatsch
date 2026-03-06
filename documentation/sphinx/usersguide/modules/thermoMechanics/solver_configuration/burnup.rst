Burnup
======

The ``burnup`` model handles the evolution of fuel burnup in the
thermo-mechanics module (OFFBEAT / GeN-Foam solid region).

The burnup field is named ``Bu`` and is expressed in
``MWd/MT_oxide``.

-------------------------------------------------------------------------------

.. rubric:: Selection in ``solverDict``

The burnup model is selected in ``solverDict`` via the ``burnup`` dictionary
(located in ``constant/`` for OFFBEAT, or in ``constant/<regionName>/`` for
GeN-Foam):

.. code-block:: cpp

   burnup
   {
       type Lassmann;
       // additional options (if any)
   }

If the ``burnup`` dictionary is not present in ``solverDict``, the burnup
evolution is neglected and the burnup field is not created.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following burnup models:

.. toctree::
   :maxdepth: 1
   :caption: Burnup models
   :glob:

   ../../../../../cppapi/generated/offbeatLib/burnup/*