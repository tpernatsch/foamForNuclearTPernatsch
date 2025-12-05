.. _post_processing:

==============
Postprocessing
==============

Postprocessing can be performed using ``paraFoam``, the standard post-processing
tool used with OpenFOAM. ``paraFoam`` is launched using the command line:

.. code :: bash

    paraFoam -region regionName


where the ``regionName`` is e.g ``fluidRegion`` or ``neutroRegion``.

Please notice that ``paraFoam`` is essentially an extension of ParaView and it
requires having ParaView installed. In the `openfoam.com <https://openfoam.com>`_
distribution, ParaView is not distributed with OpenFOAM, but needs to be
installed separately (see the `ParaView website <https://www.paraview.org/>`_). In
Ubuntu, it is normally enough to type in the terminal:

.. code :: bash

    sudo apt-get -y install paraview


In case of parallel calculations, one should first reconstruct each one of the
meshes using the command

.. code :: bash

    reconstructPar -region regionName
    # or
    reconstructPar -allRegions


where the ``regionName`` is once again e.g ``fluidRegion`` or ``neutroRegion``.
Notice that the ``<timeStep>/uniform`` folder is not a region and requires the
``-allRegions`` flag to be merged.

Besides ``paraFoam``, GeN-Foam also creates, in the case folder (or in the
processor folder for parallel simulations), the ``GeN-Foam.dat`` file that
summarizes few main quantity of interests: time(s), keff(-), power(W), flux0
(m-2s-1), TFuel_Max, TFuel_Avg TFuel_Min, TCladding_Max, TCladding_Avg,
TCladding_Min.

Useful information is also stored in the log file. The log file can be created
by adding the ``| tee log.GeN-Foam`` command to the launch command, i.e.:

.. code :: bash

    GeN-Foam | tee log.GeN-Foam
    # or
    mpirun -np <nProcessors> GeN-Foam -parallel | tee log.GeN-Foam


Python can effectively be used to extract information from the ``log.GeN-Foam``
(several examples are available in the tutorials).

In addition, OpenFOAM allows to use
`function objects <https://www.openfoam.com/documentation/guides/latest/doc/guide-function-objects.html>`_
to extract specific information after or during the simulation. Also in this
case it is essential to indicate which region you want the function object to be
applied to, for instance:

.. code :: bash

    postProcess -func singleGraph -region neutroRegion


Function objects can also employed at run time via the *controlDict* (see for
instance `2D_FFTF <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/2D_FFTF/rootCase/system/controlDict>`_).

FoamForNuclear provides some additional function objects that have been found
to be particularly useful in the nuclear field, including for example 
:doc:`generated/functionObjects/massFlow/massFlow`, 
:doc:`generated/functionObjects/TBulk/TBulk` and
:doc:`generated/functionObjects/pressureDrop/pressureDrop`.
See `here <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/src/functionObjects?ref_type=heads>`
for a full list
