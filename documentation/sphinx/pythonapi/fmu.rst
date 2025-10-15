.. _pythonapi_fmi:

------------
FMI Coupling
------------

FMU Container
-------------

FMU containers are special wrappers to simplify the manipulation of FMU with
OpenFOAM cases. Three containers are available, each based on different Python
packages:

- :class:`foamForNuclear.FMPyContainer` based on FMPy
- :class:`foamForNuclear.OMSimulatorContainer` based on OMSimulator
- :class:`foamForNuclear.PyFMIContainer` based on PyFMI

N.B: `FMPyContainer` is recommended as it offers more capabilities than the
other models.

Examples of use are provided in the following examples:

- tutorials/fmuCases/2D_PKCoupleFMI


Export OpenFOAM cases in FMU
----------------------------

FoamForNuclear can be used to export OpenFOAM cases into FMU using the
``generateCaseAsFMU`` function.

.. code-block:: python

    # Export GeN-Foam as an FMU
    ffn.generateCaseAsFMU("rootCase", "reactor")

    # >>> generate FMU named "reactor.fmu"


Multiple FMU in Parallel
------------------------

FoamForNuclear offers the possibility to run multiple couple FMUs in sequential
or in parallel.

- :class:`foamForNuclear.FMPyMasterRunner` runs in sequential
- :class:`foamForNuclear.FMPyMasterRunnerParallel` runs in parallel

These objects allow to run FMUs using a semi-implicit scheme and manage time
steps based on the smallest among all FMUs.

Example of use can be found here:

- tutorials/fmuCases/2D_PKCoupleFMIas2FMUs
