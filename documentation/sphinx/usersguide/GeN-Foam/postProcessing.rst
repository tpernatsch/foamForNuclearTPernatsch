.. _postProcessing:


Visualization and postprocessing
================================

ParaFoam and ParaView
---------------------

The most common method to vizualize and postprocess simulation results in OpenFOAM is to use ``paraFoam``. ``paraFoam`` is launched using the command line:

.. code :: bash

    paraFoam

``paraFoam`` is just an extension of ParaView and it
requires having ParaView installed. In the `openfoam.com <https://openfoam.com>`_
distribution, ParaView is not distributed with OpenFOAM, but needs to be
installed separately (see the `ParaView website <https://www.paraview.org/>`_). In
Ubuntu, it is normally enough to type in the terminal:

.. code :: bash

    sudo apt-get -y install paraview


Standard paraview can also be used as:

.. code :: bash

    touch para.foam
    paraview para.foam

where ``para.foam`` can be any empty file with extension ``.foam``.

In case of parallel calculations, one should first reconstruct each one of the
meshes using the command

.. code :: bash

    reconstructPar -region regionName
    # or
    reconstructPar -allRegions


.. note::
    The ``<timeStep>/uniform`` folder is not a region and requires the ``-allRegions`` flag to be merged.


The log file
------------

Useful information is also stored in the log file. It includes details of the numerical solution (time steps, residuals, number of interations, etc), as well as some quantities of interest such as min/max temperatures in fluid solvers and multiplication factor in k-eigenvalue solvers. 
The log file can be created
by adding the ``| tee log.GeN-Foam`` command to the launch command, i.e.:

.. code :: bash

    GeN-Foam | tee log.GeN-Foam
    # or
    mpirun -np <nProcessors> GeN-Foam -parallel | tee log.GeN-Foam


Python can effectively be used to extract information from the ``log.GeN-Foam``
(several examples are available in the tutorials).

Function objects 
----------------

OpenFOAM allows to use
`function objects <https://www.openfoam.com/documentation/guides/latest/doc/guide-function-objects.html>`_
to extract specific information after a simulation. Also in this
case it is essential to indicate which region you want the function object to be
applied to, for instance:

.. code :: bash

    postProcess -func singleGraph -region neutroRegion


Function objects can also be employed at run time via the *controlDict* (see for
instance `2D_FFTF <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/2D_FFTF/rootCase/system/controlDict>`_).

foamForNuclear provides some additional function objects with respect to standard OpenFOAM that have been found
to be particularly useful in the nuclear field.

.. including for example for 
   :ref:`mass flow rates <massFlow.H>`, 
   :ref:`pressure drops <pressureDrop.H>`
   and :ref:`bulk temperatures <TBulk.H>`.
   See `here <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/src/functionObjects?ref_type=heads>`_
   for a full list.
