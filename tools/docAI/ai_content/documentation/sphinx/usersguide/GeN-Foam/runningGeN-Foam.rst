Running GeN-Foam
================

GeN-Foam is launched like any other OpenFOAM solver by executing the following
commands in a terminal (after sourcing the OpenFOAM environment):

.. code :: bash

    GeN-Foam
    # or
    mpirun -np <nProcessors> GeN-Foam -parallel

For parallel calculations, you should first decompose each region using the command:

.. code :: bash

    decomposePar -region regionName
    # or
    decomposePar -allRegions

Here, ``decomposePar`` operates according to the standard OpenFOAM
``decomposeParDict`` to be placed in ``system/(regionName)``.
Note that ``-allRegions`` should be used to copy the ``timeStep/uniform`` folder,
which is not associated with any region.

.. note ::

    All meshes must be decomposed into the same number of domains.

In the tutorials, ``Allrun`` (and sometimes ``Allrun_parallel``) bash scripts are
provided to run the case (type ``Allrun`` or ``Allrun_parallel`` in a terminal).
These scripts include mesh decomposition for parallel cases, as well as multiple
GeN-Foam runs used, for instance, to: 1) achieve a steady state; and 2) run a transient
starting from that steady state.

.. warning::
    ``decomposePar`` and ``reconstructPar`` cannot operate on FieldFields (i.e., fields of fields, such as fields of multi-node temperatures used in sub-scale power models for porous-medium simulations). This affects the parallel decomposition and reconstruction of objects such as ``T.lumpedNuclearStructure`` and ``T.nuclearFuelPin``. The main consequence is that a user cannot automatically restart a case that has been reconstructed and re-decomposed. The best option is to simply not reconstruct or decompose before a restart. If, for whatever reason, this is needed, the user may have to do it manually or through a dedicated Python script.