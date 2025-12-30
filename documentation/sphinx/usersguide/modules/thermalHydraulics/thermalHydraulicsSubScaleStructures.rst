.. _thermalHydraulicsSubScaleStructures:


Sub-scale structures
====================

Dedicated models for specific sub-scale structures are

Power models
------------

Power models, i.e., active media that can provide and subtract energy,
including:

.. list-table:: Power models
    :widths: 50 50 50

    * - :ref:`fixedPower <fixedPower>`
      - Fixed (possibly time-dependent) power
      - *1D_CHF/imposedPower*
    * - :ref:`fixedTemperature <fixedTemperature>`
      - Fixed (possibly time-dependent) temperature
      - *1D_CHF/imposedTemperature*
    * - :ref:`heatedPin <heatedPin>`
      - Heated pin, typically used for electrically heated pins
      - *2D_KNS37-L22*
    * - :ref:`nuclearFuelPin <nuclearFuelPin>`
      - Nuclear fuel pin
      - *3D_SmallESFR* and *2D_FFTF*
    * - :ref:`lumpedNuclearStructure <lumpedNuclearStructure>`
      - Lumped-parameter nuclear structure
      - *1D_thermalMSR_pointKinetics*
    * - :ref:`nuclearSteadyStatePebble <nuclearSteadyStatePebble>`
      - Steady-state model purpose made for pebble bed reactors
      - *3D_gFHR*
    * - :ref:`nuclearFuelFMU <nuclearFuelFMU>`
      - Nuclear fuel from FMU(s)
      - N/A

Heat exchanger
--------------

A heat exchanger model that is used to model the heat transfer between two
disconnected regions, for instance representing the primary and secondary
circuit (see :ref:`heatExchanger.H <heatExchanger>` and the tutorials *1D_HX* and *2D_FFTF*)


Pump
----

A pump model used to set a (possibly time-dependent) momentum source (see
:ref:`pump.H <pump>` and tutorials *2D_FFTF* and *2D_MSFR*).
