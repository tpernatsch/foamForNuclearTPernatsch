
Running GeN-Foam
================

GeN-Foam is launched like any other OpenFOAM solver, by executing the following
commands in a terminal (after sourcing the OpenFOAM environment):

.. code :: bash

    GeN-Foam
    # or
    mpirun -np <nProcessors> GeN-Foam -parallel


In case of parallel calculations, one should first decompose each one of regions using the command

.. code :: bash

    decomposePar -region regionName
    # or
    decomposePar -allRegions


where ``decomposePar`` operates according to the standard
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


.. warning::
    ``decomposePar`` and ``reconstructPar`` are not capable to work on FieldFields (i.e., fields of fields, such as fields of multi-node temperatures emplyed in sub-scale power models for porous-medium simulations). This affect the parallel decomposition and reconstruction of objects such as ``T.lumpedNuclearStructure`` and ``T.nuclearFuelPin``. The main consequence is that a user cannot automatically restart a case that has been reconstructed and re-decomposed. The best is simply not to reconstruct/decompose before a restart. If, for whatever reason, that is needed, a user may need to do that manually or through a dedicated python script.