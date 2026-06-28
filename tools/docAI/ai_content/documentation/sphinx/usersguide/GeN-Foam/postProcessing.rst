.. _postProcessing:


Visualization and postprocessing
================================

ParaFoam and ParaView
---------------------

The most common method to visualize and postprocess simulation results in OpenFOAM is to use ``paraFoam``. ``paraFoam`` is launched from the command line:

.. code :: bash

    paraFoam

``paraFoam`` is an extension of ParaView. It requires ParaView to be installed. In the `openfoam.com <https://openfoam.com>`_ distribution, ParaView is not shipped with OpenFOAM, so it must be installed separately (see the `ParaView website <https://www.paraview.org/>`_). On Ubuntu, it is normally enough to run:

.. code :: bash

    sudo apt-get -y install paraview


Standard ParaView can also be used as:

.. code :: bash

    touch para.foam
    paraview para.foam

where ``para.foam`` can be any empty file with the extension ``.foam``.

For parallel calculations, you should first reconstruct each mesh using:

.. code :: bash

    reconstructPar -region regionName
    # or
    reconstructPar -allRegions


.. note::
    The ``<timeStep>/uniform`` folder is not a region. It requires the ``-allRegions`` flag to be merged.


The log file
-------------

Useful information is also stored in the log file. It includes details of the numerical solution (time steps, residuals, number of iterations, etc.), as well as quantities of interest such as min/max temperatures in fluid solvers and the multiplication factor in k-eigenvalue solvers.

The log file can be created by adding the ``| tee log.GeN-Foam`` command to the launch command, for example:

.. code :: bash

    GeN-Foam | tee log.GeN-Foam
    # or
    mpirun -np <nProcessors> GeN-Foam -parallel | tee log.GeN-Foam


Python can be used effectively to extract information from ``log.GeN-Foam`` (several examples are available in the tutorials).

Function objects
-----------------

OpenFOAM allows you to use
`function objects <https://www.openfoam.com/documentation/guides/latest/doc/guide-function-objects.html>`_
to extract specific information after a simulation. In this case, it is also essential to specify which region the function object should be applied to, for instance:

.. code :: bash

    postProcess -func singleGraph -region neutroRegion


Function objects can also be employed at run time via the *controlDict* (see, for instance,
`2D_FFTF <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/2D_FFTF/rootCase/system/controlDict>`_).

foamForNuclear provides additional function objects compared to standard OpenFOAM. They have been found to be particularly useful in the nuclear field:

.. toctree::
   :maxdepth: 1
   :glob:

   ../../cppapi/generated/functionObjects/**

.. raw:: html

   <br><br>