.. _solutionControlGF:

================
Solution Control
================

Solution control is mostly achived using ``controlDict``, which is in GeN-Foam a slightly extended version of the one that is normally used in other OpenFOAM solvers. 

The *controlDict* dictionary
----------------------------

The ``controlDict`` file, located in the ``system/`` directory, is the central configuration dictionary for controlling the execution of any OpenFOAM-based solver, including foamForNuclear applications such as GeN-Foam. It defines global simulation parameters, time control, output settings, and optional function objects.


Purpose
~~~~~~~

The ``controlDict`` governs:

- **Simulation time control**: Start and end times, time-step size, and write intervals.
- **Output management**: How often results are written to disk and in which format.
- **Function objects**: Post-processing utilities executed during the run (see :ref:`postProcessing <postProcessing>`).


Typical Entries
~~~~~~~~~~~~~~~

A standard ``controlDict`` includes the following keys:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Entry
     - Description
   * - ``application``
     - Name of the solver to run (e.g., ``genFoam``).
   * - ``startFrom``
     - How the simulation starts: ``startTime``, ``latestTime``, or ``firstTime``.
   * - ``startTime``
     - Initial simulation time.
   * - ``stopAt``
     - When to stop: ``endTime`` or ``writeNow``.
   * - ``endTime``
     - Final simulation time.
   * - ``deltaT``
     - Time-step size.
   * - ``writeControl``
     - Controls output frequency: ``timeStep``, ``runTime``, or ``adjustableRunTime``.
   * - ``writeInterval``
     - Interval for writing results.
   * - ``purgeWrite``
     - Number of old time directories to keep.
   * - ``functions``
     - Subdictionary for function objects (e.g., probes, field sampling).


Example controlDict
~~~~~~~~~~~~~~~~~~~

.. code-block:: foam

    application     genFoam;

    startFrom       startTime;
    startTime       0;
    stopAt          endTime;
    endTime         1000;

    deltaT          1;
    writeControl    timeStep;
    writeInterval   100;
    purgeWrite      0;

    functions
    {
        residualMonitor
        {
            type            residuals;
            libs            ("libutilityFunctionObjects.so");
            write           yes;
        }
    }




The *controlDict* extensions
----------------------------

Compared to a standard OpenFOAM controlDict, it includes
a few  specialized keywords including:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Entry
     - Description
   * - ``liquidFuel``
     - Whether or not the user is simulating a ``liquidFuel`` reactor.
   * - ``adjustTimeStep``
     - Whether or not to adjust the time step based on the conditions below.
   * - ``maxDeltaT``
     - Maximum temperature increment in a time step.
   * - ``maxCo``
     - Maximum Courant number for single-phase flow solvers.
   * - ``maxPowerVariation``
     - Maximum relative power increment in a time step.
   * - ``maxCoTwoPhase``
     - Maximum Courant number for two-phase flow solvers.
   * - ``marginToPhaseChange``
     - Temperature margin before triggering phase change. 

An additional option for mesh manipulation called ``removeBaffles``. This allows to select the physics that require the creation of a ghost mesh without baffles. This ghost mesh allows for better mesh-to-mesh projections between different physics. 

Example:

.. code :: cpp

    // In system/controlDict

    removeBaffles
    {
        fluidRegion     true;
    }

.. warning ::
    The parallel execution of the removeBaffles option is not fully tested. The user is encouraged to verify the correctness of the projection process. 


relatively complete examples of *controlDict* for single-phase flow can be found in
`2D_FFTF <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/2D_FFTF/rootCase/system/controlDict>`_
and
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/system/controlDict>`_,
while an explanation of the two-phase flow options can be found in
`1D_boiling <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/1D_boiling/system/controlDict>`_.

