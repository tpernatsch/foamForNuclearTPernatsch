# Subcritical Point-Kinetics Solver

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-onePhase-blue.svg)]() [![badge](https://img.shields.io/badge/Neutronics-pointKinetics-blue.svg)]() [![badge](https://img.shields.io/badge/Multiphysics-looseCoupling-orange.svg)]()

Author: 2022/03/28, Thomas Guilbaud, EPFL/Transmutex SA

---

This solver extends the existing GeN-Foam Point-Kinetics solver by adding an external source. Its purpose is to simulate an Accelerator-Driven System with a powerful spallation source of a few MW.


## Executable and Scripts

```bash
./Allrun_rampSource
# or
./Allrun_rampReactivity
```

To compare with other methods (make sure that the physics parameters are the same in the Python script and the GeN-Foam inputs):
```bash
python3 Allplot_rampSource.py transient/log.GeN-Foam
# or
python3 Allplot_rampReactivity.py transient/log.GeN-Foam
```

To process the files for plotting:
```bash
./parseFoam.sh transient/log.GeN-Foam
```

The script outputs:
- `file.mod`: external source modulation
- `file.power`: fission, decay, and external source power
- `file.reactivity`: all reactivities
- `file.temperature`: integrated temperature of the fuel, cladding, and coolant

To plot other transient parameters:
```bash
# Power
gnuplot -e "file='transient/log.GeN-Foam.power'" plottingTools/powerPlot.r - &> /dev/null &
# Temperature
gnuplot -e "file='transient/log.GeN-Foam.temperature'" plottingTools/temperaturePlot.r - &> /dev/null &
# Reactivity
gnuplot -e "file='transient/log.GeN-Foam.reactivity'" plottingTools/reactivityPlot.r - &> /dev/null &
```

## Results

### Comparison of GeN-Foam for a source ramp variation

![](images/2D_onePhaseAndSubcriticalPointKineticsCoupling_2023-06-01_figRampSource.png)
*Fig 1: Power evolution for a source ramp variation.*

| Method              | Value [W]    | Error [%]    |
|:--------------------|:------------:|:------------:|
| Infinite analytic   | 2.000000e+07 | -            |
| Last value analytic | 1.996341e+07 | -            |
|            GeN-Foam | 1.996340e+07 | 5.054466e-05 |
|       Runge-Kutta 4 | 1.996342e+07 | 3.625729e-05 |
|  Euler Substitution | 1.996341e+07 | 4.977224e-06 |
|  Matrix Exponential | 1.996342e+07 | 3.625729e-05 |


### Comparison of GeN-Foam for a reactivity ramp insertion

![](images/2D_onePhaseAndSubcriticalPointKineticsCoupling_2023-06-01_figRampReactivity.png)
*Fig 2: Power evolution for a reactivity ramp variation.*

| Method              | Value [W]    | Error [%]    |
|:--------------------|:------------:|:------------:|
| Infinite analytic   | 2.000000e+07 | -            |
| Last value analytic | 1.992711e+07 | -            |
|            GeN-Foam | 1.992710e+07 | 4.351604e-05 |
|       Runge-Kutta 4 | 1.992711e+07 | 1.023480e-06 |
|  Euler Substitution | 1.992710e+07 | 6.764785e-05 |
|  Matrix Exponential | 1.992711e+07 | 1.023479e-06 |