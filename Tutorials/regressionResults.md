# GeN-Foam regression test suite

Regression test of 11-oct.-2024 - 17:14:03
Tutorials successfully completed: 22/23  [![regressionTest](https://img.shields.io/badge/Alltest-failed-red.svg?style=flat-square)]()

## Summary

| Tutorial | Status | Time | GeN-Foam memory usage [Gb] |
|:---------|:------:|:-:|:-:|
| 1D_boiling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:01:51 | .076 < 0.2 |
| 1D_CHF/imposedPower | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:16:03 | .075 < 0.2 |
| 1D_CHF/imposedTemperature | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:02:58 | .075 < 0.2 |
| 1D_HX/onePhase | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:06 | .080 < 0.2 |
| 1D_HX/twoPhase | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:53 | .083 < 0.2 |
| 1D_MSR_pointKinetics | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:01:26 | .082 < 0.2 |
| 1D_PSBT_SC/Phase_Ex1_12223 | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:13 | .076 < 0.2 |
| 1D_PSBT_SC/PhaseII_Ex1_01_5215 | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:01:26 | .084 < 0.2 |
| 1D_PSBT_SC/PhaseII_Ex2_04_6770 | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 01:12:53 | .081 < 0.2 |
| 1D_thermalMSR_pointKinetics | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:49 | .081 < 0.2 |
| 2D_cavityBoussinesq | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:35 | .077 < 0.2 |
| 2D_externalSourceDiffusion | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:09 | .078 < 0.2 |
| 2D_FFTF | [![Alltest](https://img.shields.io/badge/Alltest-failed-red.svg?style=flat-square)]() |  |  < 0.3 |
| 2D_fullCoupling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:09 | .079 < 0.2 |
| 2D_KNS37-L22 | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:26:24 | .086 < 0.3 |
| 2D_MSFR | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:10:49 | .176 < 0.3 |
| 2D_onePhaseAndPointKineticsCoupling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:33 | .079 < 0.2 |
| 2D_onePhaseAndSubcriticalPointKineticsCoupling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:57 | .079 < 0.2 |
| 2D_voidMotionNoPhaseChange | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:38 | .080 < 0.2 |
| 3D_gFHR | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:13:56 | .786 < 1.5 |
| 3D_NTPfuelAssembly | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 01:28:21 | .092 < 0.2 |
| 3D_SmallESFR_NewSolverVerification | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:07:07 | .174 < 1.0 |
| Godiva_SN | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:38:01 | 1.072 < 1.5 |

--------------------------------------------------------------------------------

## 1D_boiling

Run tutorial 1D_boiling ...

Done running tutorial 1D_boiling ...

1D_boiling has converged


Time = 00:01:51

Memory = .076 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 1D_CHF/imposedPower

Run tutorial 1D_CHF imposedPower ...

Done running tutorial 1D_CHF imposedPower ...

1D_CHF imposedPower has converged

Structure temperature at the top of the channel at 90 seconds with 0.01 relative error
Perfect match:
|          | Simulated | Expected |
|:---------|:---------:|:--------:|
| TwallTop | 2610.24 | 2617.31 |

Time = 00:16:03

Memory = .075 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 1D_CHF/imposedTemperature

Run tutorial 1D_CHF imposedTemperature ...

Done running tutorial 1D_CHF imposedTemperature ...

1D_CHF imposedTemperature has converged

Heat flux of the structure at the top of the channel at 40 seconds with 0.01 relative error
Perfect match:
|           | Simulated | Expected |
|:----------|:---------:|:--------:|
| Qwall Top | 32947.6 | 32947.6 |

Time = 00:02:58

Memory = .075 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 1D_HX/onePhase

Run tutorial 1D_HX onePhase ...

Done running tutorial 1D_HX onePhase ...

1D_HX onePhase has converged


Time = 00:00:06

Memory = .080 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 1D_HX/twoPhase

Run tutorial 1D_HX twoPhase ...

Done running tutorial 1D_HX twoPhase ...

1D_HX twoPhase has converged


Time = 00:00:53

Memory = .083 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 1D_MSR_pointKinetics

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

Time = 00:01:26

Memory = .082 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 1D_PSBT_SC/Phase_Ex1_12223

Run tutorial 1D_PSBT_SC ...

Done running tutorial 1D_PSBT_SC ...

1D_PSBT_SC has converged

alpha.vapour at the end of the channel with 0.05 relative error
Perfect match:
|              | Simulated | Expected |
|:-------------|:---------:|:--------:|
| alpha vapour | 0.133727 | 0.131886 |

Time = 00:00:13

Memory = .076 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 1D_PSBT_SC/PhaseII_Ex1_01_5215

Run tutorial PSBT_Rod_Bundle_CHF PhaseII_Ex1_01_5215 ...

Done running tutorial PSBT_Rod_Bundle_CHF PhaseII_Ex1_01_5215

PSBT_Rod_Bundle_CHF PhaseII_Ex1_01_5215 has converged

alpha.vapour at the end of the channel with 0.01 relative error
Perfect match:
|                       | Simulated | Expected |
|:----------------------|:---------:|:--------:|
| Max alpha vapour      | 0.123835 | 0.123835 |
| T at max alpha vapour | 620.178 | 620.178 |

Time = 00:01:26

Memory = .084 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 1D_PSBT_SC/PhaseII_Ex2_04_6770

Run tutorial PSBT_Rod_Bundle_CHF PhaseII_Ex2_04_6770 ...

Done running tutorial PSBT_Rod_Bundle_CHF PhaseII_Ex2_04_6770

PSBT_Rod_Bundle_CHF PhaseII_Ex2_04_6770 has converged

alpha.vapour at the end of the channel with 0.05 relative error
Perfect match:
|                       | Simulated | Expected |
|:----------------------|:---------:|:--------:|
| Max alpha vapour      | 0.322952 | 0.319074 |
| T at max alpha vapour | 625.839 | 625.826 |

Time = 01:12:53

Memory = .081 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 1D_thermalMSR_pointKinetics

Run tutorial 1D_MSR_pointKinetics ...

Done running tutorial 1D_MSR_pointKinetics ...

Steady state has converged

Transient has converged

Transient end has converged

Power at the end of transient with 0.01 relative error
Perfect match:
|       | Simulated | Expected |
|:------|:---------:|:--------:|
| Power | 415251000 | 415258000 |

Time = 00:00:49

Memory = .081 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 2D_cavityBoussinesq

Run tutorial 2D_cavityBoussinesq ...

Done running tutorial 2D_cavityBoussinesq ...

2D_cavityBoussinesq has converged


Time = 00:00:35

Memory = .077 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 2D_externalSourceDiffusion

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

Time = 00:00:09

Memory = .078 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 2D_FFTF


Time =

Memory =  Gb (, expected: 0.3 Gb)

--------------------------------------------------------------------------------

## 2D_fullCoupling

Run tutorial 2D_fullCoupling ...

Done running tutorial 2D_fullCoupling ...

Steady State has converged

Transient (no driveline) has converged


Time = 00:00:09

Memory = .079 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 2D_KNS37-L22

Run tutorial 2D_KNS37-L22 ...

Done running tutorial 2D_KNS37-L22 ...

2D_KNS37-L22 has converged

Test with 0.2 relative error
Perfect match:
|                   | Simulated | Expected |
|:------------------|:---------:|:--------:|
| Time              | 12.410400262564377 | 12.514 |
| alpha.vapour      | 0.166910 | 0.178709 |
| T.activeStructure | 828.177800 | 830.2673 |

Time = 00:26:24

Memory = .086 Gb (Passed, expected: 0.3 Gb)

--------------------------------------------------------------------------------

## 2D_MSFR

Run tutorial 2D_MSFR ...

Done running tutorial 2D_MSFR ...

Thermal-hydraulic steady state has converged

Energy steady state has converged

Transient has converged

Test keff and power with 0.0001 relative error
Perfect match:
|             | Simulated | Expected |
|:------------|:---------:|:--------:|
| keff        | 0.960311 | 0.960283 |
| Total power | 20000000 | 20000000 |

Time = 00:10:49

Memory = .176 Gb (Passed, expected: 0.3 Gb)

--------------------------------------------------------------------------------

## 2D_onePhaseAndPointKineticsCoupling

Run tutorial 2D_onePhaseAndPointKineticsCoupling ...

Done running tutorial 2D_onePhaseAndPointKineticsCoupling ...

Steady State has converged

Transient (no driveline) has converged

Transient (with driveline) has converged

Transient (with boron) has converged


Time = 00:00:33

Memory = .079 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 2D_onePhaseAndSubcriticalPointKineticsCoupling

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

Time = 00:00:57

Memory = .079 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 2D_voidMotionNoPhaseChange

Run tutorial 2D_voidMotionNoPhaseChange ...

Done running tutorial 2D_voidMotionNoPhaseChange ...

2D_voidMotionNoPhaseChange has converged


Time = 00:00:38

Memory = .080 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 3D_gFHR

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

Time = 00:13:56

Memory = .786 Gb (Passed, expected: 1.5 Gb)

--------------------------------------------------------------------------------

## 3D_NTPfuelAssembly

Run tutorial 3D_NTPfuelAssembly ...

Done running tutorial 3D_NTPfuelAssembly ...

Steady state neutronics & thermal-hydraulics has converged

Simulation & reference comparison
  -> 0.1 % maximum relative error for (keff, power neutro & fluid)
  -> 5 % maximum relative error for (mass flow, Tout)

|                         | Simulated | Expected  | Error [%] |
|:------------------------|:---------:|:---------:|:---------:|
| keff                    | 1.3194364 | 1.3194342 |   -.00016 | # Passed
| Total power neutro [MW] |   937.000 |       937 |         0 | # Passed
| Total power fluid  [MW] |   936.998 |       937 |    .00021 | # Passed
| Total mass flow  [kg/s] | 30.54932000 |      31.0 |      1.45 | # Passed
| Outlet temperature  [K] | 1929.9816 |      1972 |      2.13 | # Passed

All test passed


Time = 01:28:21

Memory = .092 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 3D_SmallESFR_NewSolverVerification

Run tutorial 3D_SmallESFR ...

Done running tutorial 3D_SmallESFR ...

Steady state - Legacy has converged

Steady state - New has converged

keff test with 0.001 relative error
Perfect match:
|                | Simulated | Expected |
|:---------------|:---------:|:--------:|
| keff - Legacy  | 0.936873  | 0.936827 |
| keff - New     | 0.936872  | 0.936827 |

Time = 00:07:07

Memory = .174 Gb (Passed, expected: 1.0 Gb)

--------------------------------------------------------------------------------

## Godiva_SN

Run tutorial Godiva_SN ...

Done running tutorial Godiva_SN ...

Godiva_SN has converged

keff with 0.001 relative error
Perfect match:
|      | Simulated | Expected |
|:-----|:---------:|:--------:|
| keff | 0.980106 | 0.980106 |

Time = 00:38:01

Memory = 1.072 Gb (Passed, expected: 1.5 Gb)
