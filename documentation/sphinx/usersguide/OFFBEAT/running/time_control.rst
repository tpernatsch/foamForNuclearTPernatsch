Time control
==============================

This page documents the time-step control and write-time control mechanisms available in
OFFBEAT.


User time units
---------------

OFFBEAT supports user-selected time units in addition to the default seconds used by most
OpenFOAM solvers. The time unit can be set in ``controlDict`` to ``seconds``, ``hours``, or
``days``.

.. warning::

    All time-dependent lists (e.g. LHGR histories, fast flux histories) and all time-step
    criteria expressed in time units must be provided in user time. OFFBEAT converts these
    values internally to seconds.

.. code-block:: c

    // User time units (all lists and criteria below use this unit)
    userTime
    {
        type seconds; // alternatives: hours; days;
    }


Enabling time-step control
--------------------------

Automatic time-step adjustment is enabled with:

.. code-block:: c

    adjustTimeStep yes;

.. warning::

    The legacy keyword ``adjustableTimeStep`` is still recognised for backward compatibility
    but is deprecated. If present, OFFBEAT prints a warning.

When the switch is disabled, all time-step criteria described below are ignored and the
simulation advances with the fixed ``deltaT`` specified in ``controlDict``.


Purpose of time-step control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

During transients, abrupt changes in power, burnup, fission gas release, porosity transport,
creep, or plasticity can lead to instability if the time step is too large. The time
controller monitors selected quantities and proposes a new time step when one or more
configured thresholds are exceeded. The proposed value is then constrained by absolute and
relative limits (e.g. ``minDeltaT``, ``maxDeltaT``, ``maxRelativeDeltaTIncrease``).

This mechanism results in smaller time steps during rapid physical changes and larger time
steps when the evolution is smooth.

.. note::

    Not all keywords listed below are required for every simulation. Criteria associated
    with a specific physics are considered only when that physics is active.

.. warning::

    Stricter thresholds typically increase the number of time steps and can significantly
    raise computational cost.


Available options
-----------------

All entries in this section are considered only when ``adjustTimeStep yes;`` is set.
Defaults are indicated where applicable. The value ``GREAT`` denotes “no practical limit”.


General :math:`\Delta t` limits
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``maxDeltaT``  
  Maximum allowed time step in user time units. Default: ``GREAT``.

* ``minDeltaT``  
  Minimum allowed time step in user time units. Default: ``1`` (in user time units).

* ``maxRelativeDeltaTIncrease``  
  Maximum factor by which :math:`\Delta t` may increase from one step to the next.
  Default: ``GREAT``.


Burnup criterion
~~~~~~~~~~~~~~~~

* ``maxBurnupIncrease``  
  Maximum allowed burnup increment per step in MWd/kgU. Internally converted to MWd/t.
  Default: ``GREAT``.


Fission gas release criteria
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Legacy keyword (deprecated):

* ``maxFGR``  
  Deprecated limit on FGR growth, expressed as total moles per time step.

Current interface:

* ``maxLocalFgrChange``  
  Maximum allowed local FGR increment per step, in percent (not as a fraction).
  Default: ``GREAT``.

* ``maxTotalFgrChange``  
  Maximum allowed total FGR increment per step, in percent (not as a fraction).
  Default: ``GREAT``.

If ``maxFGR`` is present, OFFBEAT uses the legacy behaviour for backward compatibility. This
behaviour is deprecated and should be migrated to the current interface.


Power-based criteria
~~~~~~~~~~~~~~~~~~~~

* ``maxRelativePowerIncrease``  
  Maximum allowed relative increase of total power per step. Default: ``GREAT``.

* ``maxRelativePowerDecrease``  
  Maximum allowed relative decrease of total power per step. Default: ``GREAT``.

* ``maxAbsolutePowerIncrease``  
  Maximum allowed absolute increase of power density (W/m3). For fuel rods,
  ``maxAbsoluteLhgrIncrease`` is generally preferred, since it is difficult to define a
  robust power-density criterion. Default: ``GREAT``.


LHGR criterion
~~~~~~~~~~~~~~

* ``maxAbsoluteLhgrIncrease``  
  Maximum permitted change in linear heat generation rate (W/m) per step. Default: ``GREAT``.


Porosity-transport Courant limit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``maxCourantNoPorosity``  
  Maximum allowed Courant number associated with porosity transport. Default: ``GREAT``.


Material-property stability criterion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``maxDeltaMatTol``  
  Maximum tolerated relative variation of material properties (e.g. :math:`\rho`, :math:`C_p`,
  :math:`k`, emissivity, :math:`E`, :math:`\nu`, :math:`\alpha`). Default: ``GREAT``.


Creep and plasticity criteria
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``maxAverageCreep``  
  Maximum allowed average creep increment per step. Default: ``GREAT``.

* ``maxMaximumCreep``  
  Maximum allowed pointwise creep increment per step. Default: ``GREAT``.

* ``maxAveragePlasticIncrement``  
  Maximum allowed average plastic strain increment per step. Default: ``GREAT``.

* ``maxMaximumPlasticIncrement``  
  Maximum allowed pointwise plastic increment per step. Default: ``GREAT``.


Write-time control
------------------

Write-time control operates independently of the time-step controller.

The preferred keyword is:

.. code-block:: c

    adjustWriteTimeStep yes;

.. warning::

    The deprecated keyword ``adjustableWriteTimeStep`` is still recognised but triggers a
    warning.

When enabled, OFFBEAT expects a list of output times:

.. code-block:: c

    writeTimeStepList ( t1 t2 t3 ... );

Entries are provided in user time units and converted internally to seconds. The solver
temporarily overrides the write interval so that output is written exactly (and only) at the
listed times.


Example ``controlDict`` snippet
-------------------------------

.. code-block:: c

    // ---------------------------------------------
    // Time-step control
    // ---------------------------------------------

    // User time units (NOTE: all time lists and time criteria
    // need to be provided in the selected time unit)
    userTime
    {
        type seconds; // hours; days;
    }

    // Main time step settings
    adjustTimeStep              yes;
    maxDeltaT                   1e5;  // user time units
    minDeltaT                   1;    // user time units
    maxRelativeDeltaTIncrease   1.5;

    // Burnup limits in MWd/kgU
    maxBurnupIncrease           0.5;

    // FGR limits in % (NOT as a fraction)
    maxLocalFgrChange           1;
    maxTotalFgrChange           5;

    // Power behaviour
    maxRelativePowerIncrease    0.05;
    maxRelativePowerDecrease    0.10;
    // For fuel rods, maxAbsoluteLhgrIncrease is preferred
    // maxAbsolutePowerIncrease  1e7;

    // LHGR behaviour in W/m
    maxAbsoluteLhgrIncrease     100;

    // Porosity transport
    maxCourantNoPorosity        0.8;

    // Material stability
    maxDeltaMatTol              0.05;

    // Creep and plasticity
    maxAverageCreep             1e-4;
    maxMaximumCreep             5e-4;
    maxAveragePlasticIncrement  1e-4;
    maxMaximumPlasticIncrement  5e-4;

    // ---------------------------------------------
    // Write-time control
    // ---------------------------------------------
    adjustWriteTimeStep         yes;
    writeTimeStepList           (1 20 100);


Deprecated keywords
-------------------

Deprecated keywords remain supported for legacy cases but should be migrated to the
preferred forms.

* ``adjustableTimeStep``  
  Supported, deprecated. Replacement: ``adjustTimeStep``.

* ``maxFGR``  
  Supported, deprecated. Replacement: ``maxLocalFgrChange`` and ``maxTotalFgrChange``.

* ``adjustableWriteTimeStep``  
  Supported, deprecated. Replacement: ``adjustWriteTimeStep``.
