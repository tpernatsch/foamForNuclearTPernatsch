# 1D Critical Heat Flux

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-twoPhase-blue.svg)]()

## Description

This tutorial provides two examples of 1-D channels with boiling water and
achievement of critical heat flux (CHF) conditions. In both examples, CHF is
studied under (i) imposed power and (ii) imposed temperature.

For the imposed-power case, the power is gradually increased and then
decreased to reproduce the well-known hysteresis caused by the boiling
crisis.

For the imposed-temperature case, the temperature is increased only with
time. This demonstrates that, with imposed temperature, the simulation
passes through the Leidenfrost temperature.

The resulting heat flux versus wall temperature curves can be plotted using
the following scripts, which are available in each case folder:
``` bash
python3 printResults.py
```

**N.B.1**: The CHF and post-CHF models have not been validated yet. Use them
with care.

**N.B.2**: For the moment, the critical heat flux can only be imposed at a
constant value. No models or lookup tables have been implemented yet for its
prediction.