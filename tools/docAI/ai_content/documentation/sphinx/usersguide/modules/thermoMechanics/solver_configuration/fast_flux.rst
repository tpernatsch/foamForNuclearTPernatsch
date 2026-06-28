Fast flux
=========

The ``fastFlux`` model defines the fast neutron flux and governs the
evolution of fast fluence in the thermo-mechanics module
(OFFBEAT / GeN-Foam solid region).

Two fields are used:

* ``fastFlux`` – fast neutron flux (``n/cm^2/s``)
* ``fastFluence`` – fast neutron fluence (``n/cm^2``)

-------------------------------------------------------------------------------

.. rubric:: Selection in ``solverDict``

The fast flux model is selected in ``solverDict`` via the ``fastFlux``
dictionary (located in ``constant/`` for OFFBEAT, or in
``constant/<regionName>/`` for GeN-Foam):

.. code-block:: cpp

   fastFlux
   {
       type constant;
       // additional options (if any)
   }

If the ``fastFlux`` dictionary is not present in ``solverDict``, fast flux
and fast fluence are neglected, and the corresponding fields are not
created.

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following fast flux
models:

.. toctree::
   :maxdepth: 1
   :caption: Fast flux models
   :glob:

   ../../../../../cppapi/generated/offbeatLib/fastFlux/*