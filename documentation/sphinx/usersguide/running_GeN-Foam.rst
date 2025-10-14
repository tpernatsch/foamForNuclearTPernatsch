================
Running GeN-Foam
================

GeN-Foam is launched like any OpenFOAM solver, by executing the following
commands in a terminal (after sourcing the OpenFOAM environment):

.. code :: bash

    GeN-Foam
    # or
    mpirun -np <nProcessors> GeN-Foam -parallel


In case of parallel calculations, one should decompose each one of three meshes
using the command

.. code :: bash

    decomposePar -region regionName
    # or
    decomposePar -allRegions


where the ``regionName`` can be ``fluidRegion``, ``neutroRegion`` or
``thermalMechanicalRegion``. ``decomposePar`` operates according to the standard
OpenFOAM ``decomposeParDict`` to be placed in ``system/(regionName)``. Note that
the `-allRegions` should be used to copy the ``timeStep/uniform`` folder that is
not a region.

.. note ::

    All meshes must be decomposed in the same number of domains.

In the tutorials, ``Allrun`` (and sometimes ``Allrun_parallel``) bash scripts are
provided that can be used to run the case (one should just type ``Allrun`` or
``Allrun_parallel`` in a terminal), including mesh decomposition for parallel
cases, as well as multiple GeN-Foam runs that are used for instance to: 1)
achieve a steady-state; 2) run a transient starting from that steady-state.

More details on post-processing :ref:`post_processing <post_processing>`.
