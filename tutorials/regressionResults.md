# GeN-Foam regression test suite

Regression test of 29-Jun-2026 - 19:02:28  
OpenFOAM Version: v2606  
Tutorials successfully completed: 30/30  [![regressionTest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]()  

## Summary

| Tutorial | Status | Time | GeN-Foam memory usage [Gb] |
|:---------|:------:|:-:|:-:|
| ../tests | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:01 | 0 < 0.2 |
| guidedCases/1_reactorSlab_1D_1Gr_neutronicDiffusion | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:01 | .026 < 0.2 |
| guidedCases/2_reactorSlabReflected_1D_1Gr_neutronicDiffusion | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:02 | .059 < 0.2 |
| featureCases/1D_boiling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:21:36 | .084 < 0.2 |
| featureCases/1D_CHF/imposedPower | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:38:24 | .085 < 0.2 |
| featureCases/1D_CHF/imposedTemperature | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:07:08 | .084 < 0.2 |
| featureCases/1D_HX/onePhase | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:12 | .091 < 0.2 |
| featureCases/1D_HX/twoPhase | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:03:09 | .092 < 0.2 |
| featureCases/1D_PSBT_SC/Phase_Ex1_12223 | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:19 | .085 < 0.2 |
| featureCases/1D_PSBT_SC/PhaseII_Ex1_01_5215 | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:04:14 | .095 < 0.2 |
| featureCases/1D_PSBT_SC/PhaseII_Ex2_04_6770 | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 02:20:08 | .090 < 0.2 |
| featureCases/2D_cavityBoussinesq | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:01:25 | .087 < 0.2 |
| featureCases/2D_externalSourceDiffusion | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:03 | .087 < 0.2 |
| featureCases/2D_flowOverHeatedPlate | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:15:00 | .131 < 0.2 |
| featureCases/2D_fullCoupling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:13 | .090 < 0.2 |
| featureCases/2D_KNS37-L22 | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 04:59:22 | .098 < 0.3 |
| featureCases/2D_onePhaseAndPointKineticsCoupling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:01:16 | .090 < 0.2 |
| featureCases/2D_onePhaseAndSubcriticalPointKineticsCoupling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:02:18 | .090 < 0.2 |
| featureCases/2D_voidMotionNoPhaseChange | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:01:33 | .090 < 0.2 |
| reactorCases/1D_MSR_pointKinetics | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:03:59 | .091 < 0.2 |
| reactorCases/1D_thermalMSR_pointKinetics | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:03:50 | .091 < 0.2 |
| reactorCases/2D_FFTF | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 03:56:30 | .159 < 0.3 |
| reactorCases/2D_MSFR | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:37:29 | .207 < 0.3 |
| reactorCases/3D_genericParticleBedNTP | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 01:18:14 | .186 < 1.5 |
| reactorCases/3D_gFHR | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 01:01:01 | .990 < 1.5 |
| reactorCases/3D_HTR-10 | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 03:38:39 | 1.055 < 1.5 |
| reactorCases/3D_KIWI-B-4E | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 03:41:57 | .667 < 1.5 |
| reactorCases/3D_NTPfuelAssembly | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 01:58:20 | .138 < 0.2 |
| reactorCases/3D_SmallESFR | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:09:02 | .182 < 1.0 |
| reactorCases/Godiva_SN | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 02:38:40 | 1.336 < 1.5 |

--------------------------------------------------------------------------------

## ../tests

Test mole fraction x_H with 5% threshold ... passed with 5 points

Test rho (density) with 15% threshold ... passed with 112 points

Test cp (isobaric specific heat capacity) with 10% threshold ... passed with 40 points

Test cv (isochoric specific heat capacity) with 15% threshold ... passed with 28 points

Test mu (dynamic viscosity) with 25% threshold ... 
     WARNING: mu(1.03421e+07 Pa, 83.3333 K) = 6.65438e-06 Pa.s (exp: 5.1228e-06) -> rel error = 0.260092
     WARNING: mu(3.10264e+07 Pa, 83.3333 K) = 1.17613e-05 Pa.s (exp: 8.70808e-06) -> rel error = 0.298321
passed with 112 points

Test Ha (specific enthalpy) with 50% threshold ... passed with 112 points

Test D (diffusion) with 5% threshold ... passed with 6 points

Test S (specific entropy) with 15% threshold ... passed with 112 points

Test 1D polyharmonic spline RBF with 0.01% threshold ... passed with 4 points

Test 2D polyharmonic spline RBF with 0.01% threshold ... passed with 5 points

Test 2D polyharmonic spline RBF for small values with 0.02% threshold ... passed with 2 points


Time = 00:00:01

Memory = 0 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## guidedCases/1_reactorSlab_1D_1Gr_neutronicDiffusion

Run tutorial 1_reactorSlab_1D_1Gr_neutronicDiffusion ...

Done running tutorial 1_reactorSlab_1D_1Gr_neutronicDiffusion ...

1_reactorSlab_1D_1Gr_neutronicDiffusion has converged 

keff with 1e-05 relative error
|      | Simulated | Expected |
|:-----|:---------:|:--------:|
| keff | 1.0000037 |      1.0 |

All tests passed

Time = 00:00:01

Memory = .026 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## guidedCases/2_reactorSlabReflected_1D_1Gr_neutronicDiffusion

Run tutorial 2_reactorSlabReflected_1D_1Gr_neutronicDiffusion ...

Done running tutorial 2_reactorSlabReflected_1D_1Gr_neutronicDiffusion ...

2_reactorSlabReflected_1D_1Gr_neutronicDiffusion has converged 

keff with 1e-05 relative error
|      | Simulated | Expected |
|:-----|:---------:|:--------:|
| keff | 1.0000227 |      1.0 |

All tests passed

Time = 00:00:02

Memory = .059 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/1D_boiling

Run tutorial 1D_boiling ...

Done running tutorial 1D_boiling ...

1D_boiling has converged 


Time = 00:21:36

Memory = .084 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/1D_CHF/imposedPower

Run tutorial 1D_CHF imposedPower ...

Done running tutorial 1D_CHF imposedPower ...

1D_CHF imposedPower has converged 

Structure temperature at the top of the channel at 90 seconds with 0.01 relative error
Perfect match:
|          | Simulated | Expected |
|:---------|:---------:|:--------:|
| TwallTop | 2611.15 | 2617.31 |

Time = 00:38:24

Memory = .085 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/1D_CHF/imposedTemperature

Run tutorial 1D_CHF imposedTemperature ...

Done running tutorial 1D_CHF imposedTemperature ...

1D_CHF imposedTemperature has converged 

Heat flux of the structure at the top of the channel at 40 seconds with 0.01 relative error
Perfect match:
|           | Simulated | Expected |
|:----------|:---------:|:--------:|
| Qwall Top | 32984.2 | 32947.6 |

Time = 00:07:08

Memory = .084 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/1D_HX/onePhase

Run tutorial 1D_HX onePhase ...

Done running tutorial 1D_HX onePhase ...

1D_HX onePhase has converged 


Time = 00:00:12

Memory = .091 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/1D_HX/twoPhase

Run tutorial 1D_HX twoPhase ...

Done running tutorial 1D_HX twoPhase ...

1D_HX twoPhase has converged 


Time = 00:03:09

Memory = .092 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/1D_PSBT_SC/Phase_Ex1_12223

Run tutorial 1D_PSBT_SC ...

Done running tutorial 1D_PSBT_SC ...

1D_PSBT_SC has converged 

alpha.vapour at the end of the channel with 0.05 relative error
Perfect match:
|              | Simulated | Expected |
|:-------------|:---------:|:--------:|
| alpha vapour | 0.133727 | 0.131886 |

Time = 00:00:19

Memory = .085 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/1D_PSBT_SC/PhaseII_Ex1_01_5215

Run tutorial PSBT_Rod_Bundle_CHF PhaseII_Ex1_01_5215 ...

Done running tutorial PSBT_Rod_Bundle_CHF PhaseII_Ex1_01_5215 

PSBT_Rod_Bundle_CHF PhaseII_Ex1_01_5215 has converged 

alpha.vapour at the end of the channel with 0.01 relative error
Perfect match:
|                       | Simulated | Expected |
|:----------------------|:---------:|:--------:|
| Max alpha vapour      | 0.123835 | 0.123835 |
| T at max alpha vapour | 620.178 | 620.178 |

Time = 00:04:14

Memory = .095 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/1D_PSBT_SC/PhaseII_Ex2_04_6770

Run tutorial PSBT_Rod_Bundle_CHF PhaseII_Ex2_04_6770 ...

Done running tutorial PSBT_Rod_Bundle_CHF PhaseII_Ex2_04_6770 

PSBT_Rod_Bundle_CHF PhaseII_Ex2_04_6770 has converged 

alpha.vapour at the end of the channel with 0.05 relative error
Perfect match:
|                       | Simulated | Expected |
|:----------------------|:---------:|:--------:|
| Max alpha vapour      | 0.322955 | 0.319074 |
| T at max alpha vapour | 625.839 | 625.826 |

Time = 02:20:08

Memory = .090 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/2D_cavityBoussinesq

Run tutorial 2D_cavityBoussinesq ...

Done running tutorial 2D_cavityBoussinesq ...

2D_cavityBoussinesq has converged 


Time = 00:01:25

Memory = .087 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/2D_externalSourceDiffusion

Run tutorial 2D_externalSourceDiffusion ...

Done running tutorial 2D_externalSourceDiffusion ...

Steady-state has converged 

The relation S/(1-keff) = S + nuSigmaF*phi has been verified

Test power and fluxes  
Perfect match:
|             | Simulated | Expected |
|:------------|:---------:|:--------:|
| Total power | 53214.3 | 53214.3 |
| Flux 0      | 23260.9 | 23260.9 |
| Flux 1      | 29953.4 | 29953.4 |
| Source 0    | 3346.28 | 3346.28 |
| Source 1    | 3346.28 | 3346.28 |

Time = 00:00:03

Memory = .087 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/2D_flowOverHeatedPlate

Running tutorials...
Done running tutorials.

Results comparison:
| Case              | Tmax (sim) | Tmax (exp) | Converged | Status |
|:------------------|:----------:|:----------:|:---------:|:------:|
| CHTLoop           | 309.278000 | 309.287000 |       YES |     OK |
| boundaryCoupling  | 309.278000 | 309.287000 |       YES |     OK |

Time = 00:15:00

Memory = .131 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/2D_fullCoupling

Run tutorial 2D_fullCoupling ...

Done running tutorial 2D_fullCoupling ...

Steady State has converged 

Transient (no driveline) has converged 


Time = 00:00:13

Memory = .090 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/2D_KNS37-L22

Run tutorial 2D_KNS37-L22 ...

Done running tutorial 2D_KNS37-L22 ...

2D_KNS37-L22 has converged 

Test with 0.2 relative error
Perfect match:
|                   | Simulated | Expected |
|:------------------|:---------:|:--------:|
| Time              | 12.484178273178323 | 12.514 |
| alpha.vapour      | 0.169341 | 0.178709 |
| T.activeStructure | 828.801800 | 830.2673 |

Time = 04:59:22

Memory = .098 Gb (Passed, expected: 0.3 Gb)

--------------------------------------------------------------------------------

## featureCases/2D_onePhaseAndPointKineticsCoupling

Run tutorial 2D_onePhaseAndPointKineticsCoupling ...

Done running tutorial 2D_onePhaseAndPointKineticsCoupling ...

Steady State has converged 

Transient (no driveline) has converged 

Transient (with driveline) has converged 

Transient (with boron) has converged 


Time = 00:01:16

Memory = .090 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/2D_onePhaseAndSubcriticalPointKineticsCoupling

Run tutorial 2D_onePhaseAndSubcriticalPointKineticsCoupling rampSource ...

Run tutorial 2D_onePhaseAndSubcriticalPointKineticsCoupling rampReactivity ...

Done running tutorial 2D_onePhaseAndSubcriticalPointKineticsCoupling ...

Steady-State has converged 

Transient ramp source has converged 

Transient ramp reactivity has converged 

Final powers with 0.001 relative error

Perfect match:
|                             | Simulated | Expected | Theory |
|:----------------------------|:---------:|:--------:|:------:|
| Total power ramp Source     | 19963400.000000 | 19963400 | 20000000 |
| Total power ramp Reactivity | 19927100.000000 | 19927100 | 20000000 |

Time = 00:02:18

Memory = .090 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/2D_voidMotionNoPhaseChange

Run tutorial 2D_voidMotionNoPhaseChange ...

Done running tutorial 2D_voidMotionNoPhaseChange ...

2D_voidMotionNoPhaseChange has converged 


Time = 00:01:33

Memory = .090 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## reactorCases/1D_MSR_pointKinetics

Run tutorial 1D_MSR_pointKinetics ...

Done running tutorial 1D_MSR_pointKinetics ...

Steady state has converged 

Transient has converged 

Transient end has converged 

Power at the end of transient with 0.01 relative error
Perfect match:
|       | Simulated | Expected |
|:------|:---------:|:--------:|
| Power | 613807000 | 613794000 |

Time = 00:03:59

Memory = .091 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## reactorCases/1D_thermalMSR_pointKinetics

Run tutorial 1D_MSR_pointKinetics ...

Done running tutorial 1D_MSR_pointKinetics ...

Steady state has converged 

Transient has converged 

Transient end has converged 

Transient substepping has converged 

Transient end substepping has converged 

Power at the end of transient with 0.01 relative error
Perfect match:
|       | Simulated | Expected |
|:------|:---------:|:--------:|
| Power | 415251000 | 415258000 |
Power at the end of transient_substepping with 0.01 relative error
Perfect match:
|       | Simulated | Expected |
|:------|:---------:|:--------:|
| Power | 415249000 | 415258000 |

Time = 00:03:50

Memory = .091 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## reactorCases/2D_FFTF

Run tutorial 2D_FFTF ...

Done running tutorial 2D_FFTF ...

Energy steady state has converged 

Transient has converged 

Reactivity contributions with 0.01 relative error
Perfect match:
|             | Simulated | Expected |
|:------------|:---------:|:--------:|
| Doppler     | 238.153 | 238.1535 |
| Cladding    | 0.8704383 | 0.8706593 |
| Density     | -2.522382 | -2.523688 |
| Structures  | -5.553116 | -5.55205 |
| Driveline   | -7.249734 | -7.249654 |
| GEM         | -464.8492 | -464.8492 |
| Total power | 28762.06 | 28761.75 |

Time = 03:56:30

Memory = .159 Gb (Passed, expected: 0.3 Gb)

--------------------------------------------------------------------------------

## reactorCases/2D_MSFR

Run tutorial 2D_MSFR ...

Done running tutorial 2D_MSFR ...

Thermal-hydraulic steady state has converged 

Energy steady state has converged 

Transient has converged 

Test keff and power with 0.001 relative error  
Perfect match:
|             | Simulated | Expected |
|:------------|:---------:|:--------:|
| keff        | 0.959876 | 0.960283 |
| Total power | 20000000 | 20000000 |

Time = 00:37:29

Memory = .207 Gb (Passed, expected: 0.3 Gb)

--------------------------------------------------------------------------------

## reactorCases/3D_genericParticleBedNTP

Run tutorial 3D_genericParticleBedNTP ...

Done running tutorial 3D_genericParticleBedNTP ...

3D_genericParticleBedNTP has converged 

|                 Quantity                |   Simulated  | Expected  |
|:----------------------------------------|:------------:|:---------:|
| innerPropellantOutlet15                 | 805.29712    | 804.86569 |
| innerPropellantOutlet16                 | 1362.5251    | 1361.7657 |
| innerPropellantOutlet17                 | 1580.1392    | 1579.3129 |
| innerPropellantOutlet18                 | 1577.9273    | 1577.1483 |
| innerPropellantOutlet22                 | 1151.9201    | 1151.2722 |
| innerPropellantOutlet23                 | 1435.9134    | 1435.1303 |

Perfect match

Time = 01:18:14

Memory = .186 Gb (Passed, expected: 1.5 Gb)

--------------------------------------------------------------------------------

## reactorCases/3D_gFHR

Run tutorial 3D_gFHR ...

Done running tutorial 3D_gFHR ...

Steady state NSSP has converged 

Steady state LPS has converged 

Max temperatures (avg min max) in trisos  
Perfect match:
|           | Simulated | Expected |
|:----------|:---------:|:--------:|
| Tfmax avg | 981.812 | 981.808 |
| Tfmax min | 900.664 | 900.664 |
| Tfmax max | 1062.96 | 1062.87 |

Time = 01:01:01

Memory = .990 Gb (Passed, expected: 1.5 Gb)

--------------------------------------------------------------------------------

## reactorCases/3D_HTR-10

Run tutorial 3D_HTR-10 ...

Done running tutorial 3D_HTR-10 ...

Steady state has converged 

Max temperatures (avg min max) in trisos
Perfect match:
|           | Simulated | Expected |
|:----------|:---------:|:--------:|
| Tfmax avg | 774.452 | 774.446 |
| Tfmax min | 523.15 | 523.15 |
| Tfmax max | 1441.38 | 1441.35 |

Time = 03:38:39

Memory = 1.055 Gb (Passed, expected: 1.5 Gb)

--------------------------------------------------------------------------------

## reactorCases/3D_KIWI-B-4E

Run tutorial 3D_KIWI-B-4E ...

Done running tutorial 3D_KIWI-B-4E ...

Neutronics steady state has converged 

Fluid-Mechanics steady state has converged 

Thermal-hydraulic steady state has converged 

Simulation & reference comparison
  -> 0.001 % maximum relative error for (keff, power neutro & fluid)
  -> 0.05 % maximum relative error for (Tout)

|                         | Simulated | Expected  | Error [%] |
|:------------------------|:---------:|:---------:|:---------:|
| keff                    | 1.0940776 | 1.0940776 |         0 | # Passed
| Total power neutro [MW] |     905.0 |       905 |         0 | # Passed
| Total power fluid  [MW] | 905.00011 |       905 |   -.00001 | # Passed
| Outlet temperature  [K] |  1762.667 |  1762.667 |         0 | # Passed

All test passed


Time = 03:41:57

Memory = .667 Gb (Passed, expected: 1.5 Gb)

--------------------------------------------------------------------------------

## reactorCases/3D_NTPfuelAssembly

Run tutorial 3D_NTPfuelAssembly ...

Done running tutorial 3D_NTPfuelAssembly ...

Steady state neutronics & thermal-hydraulics has converged 

Simulation & reference comparison
  -> 0.1 % maximum relative error for (keff, power neutro & fluid)
  -> 5 % maximum relative error for (mass flow, Tout)

|                         | Simulated | Expected  | Error [%] |
|:------------------------|:---------:|:---------:|:---------:|
| keff                    | 1.3194348 | 1.3194342 |   -.00004 | # Passed
| Total power neutro [MW] |   937.000 |       937 |         0 | # Passed
| Total power fluid  [MW] |   937.003 |       937 |   -.00032 | # Passed
| Total mass flow  [kg/s] | 30.54934250 |      31.0 |      1.45 | # Passed
| Outlet temperature  [K] | 1929.9829 |      1972 |      2.13 | # Passed

All test passed


Time = 01:58:20

Memory = .138 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## reactorCases/3D_SmallESFR

Run tutorial 3D_SmallESFR ...

Done running tutorial 3D_SmallESFR ...

Steady state has converged 

keff test with 0.001 relative error
Perfect match:
|       | Simulated | Expected |
|:------|:---------:|:--------:|
| keff  | 0.936873  | 0.936827 |

Time = 00:09:02

Memory = .182 Gb (Passed, expected: 1.0 Gb)

--------------------------------------------------------------------------------

## reactorCases/Godiva_SN

Run tutorial Godiva_SN ...

Done running tutorial Godiva_SN ...

Godiva_SN has converged 

keff with 0.001 relative error
Perfect match:
|      | Simulated | Expected |
|:-----|:---------:|:--------:|
| keff | 0.980106 | 0.980106 |

Time = 02:38:40

Memory = 1.336 Gb (Passed, expected: 1.5 Gb)
