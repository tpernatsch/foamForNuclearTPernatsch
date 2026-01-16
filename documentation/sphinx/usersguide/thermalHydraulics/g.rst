.. _userguide_thermalhydraulics_g:

The *g* dictionary
------------------

The *g* dictionary can be found under *constant/(nameOfFluidRegion)/*. It is a standard
OpenFOAM dictionary that allows specifying the gravitational acceleration.

Example:

.. code :: cpp

   dimensions      [0 1 -2 0 0 0 0];
   value          (0 0 -9.81);
