# GeN-Foam regression test suite

Regression test of 01-nov.-2024 - 11:04:30  
Tutorials successfully completed: 23/23  [![regressionTest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]()  

## Summary

| Tutorial | Status | Time | GeN-Foam memory usage [Gb] |
|:---------|:------:|:-:|:-:|
| 1D_boiling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:01:59 | .078 < 0.2 |
| 1D_CHF/imposedPower | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:28:12 | .078 < 0.2 |
| 1D_CHF/imposedTemperature | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:03:57 | .078 < 0.2 |
| 1D_HX/onePhase | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:08 | .084 < 0.2 |
| 1D_HX/twoPhase | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:01:18 | .085 < 0.2 |
| 1D_MSR_pointKinetics | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:01:31 | .084 < 0.2 |
| 1D_PSBT_SC/Phase_Ex1_12223 | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:12 | .078 < 0.2 |
| 1D_PSBT_SC/PhaseII_Ex1_01_5215 | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:01:31 | .087 < 0.2 |
| 1D_PSBT_SC/PhaseII_Ex2_04_6770 | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 02:07:15 | .083 < 0.2 |
| 1D_thermalMSR_pointKinetics | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:59 | .084 < 0.2 |
| 2D_cavityBoussinesq | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:46 | .080 < 0.2 |
| 2D_externalSourceDiffusion | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:10 | .080 < 0.2 |
| 2D_FFTF | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 02:44:26 | .152 < 0.3 |
| 2D_fullCoupling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:17 | .082 < 0.2 |
| 2D_KNS37-L22 | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:52:05 | .088 < 0.3 |
| 2D_MSFR | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:20:47 | .179 < 0.3 |
| 2D_onePhaseAndPointKineticsCoupling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:52 | .082 < 0.2 |
| 2D_onePhaseAndSubcriticalPointKineticsCoupling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:01:33 | .082 < 0.2 |
| 2D_voidMotionNoPhaseChange | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:01:18 | .083 < 0.2 |
| 3D_gFHR | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:27:24 | .788 < 1.5 |
| 3D_NTPfuelAssembly | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 02:22:16 | .094 < 0.2 |
| 3D_SmallESFR_NewSolverVerification | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:11:20 | .175 < 1.0 |
| Godiva_SN | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 01:27:03 | 1.076 < 1.5 |

--------------------------------------------------------------------------------

## 1D_boiling

Run tutorial 1D_boiling ...

Done running tutorial 1D_boiling ...

1D_boiling has converged 


Time = 00:01:59

Memory = .078 Gb (Passed, expected: 0.2 Gb)

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

Time = 00:28:12

Memory = .078 Gb (Passed, expected: 0.2 Gb)

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

Time = 00:03:57

Memory = .078 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 1D_HX/onePhase

Run tutorial 1D_HX onePhase ...

Done running tutorial 1D_HX onePhase ...

1D_HX onePhase has converged 


Time = 00:00:08

Memory = .084 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 1D_HX/twoPhase

Run tutorial 1D_HX twoPhase ...

Done running tutorial 1D_HX twoPhase ...

1D_HX twoPhase has converged 


Time = 00:01:18

Memory = .085 Gb (Passed, expected: 0.2 Gb)

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

Time = 00:01:31

Memory = .084 Gb (Passed, expected: 0.2 Gb)

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

Time = 00:00:12

Memory = .078 Gb (Passed, expected: 0.2 Gb)

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

Time = 00:01:31

Memory = .087 Gb (Passed, expected: 0.2 Gb)

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

Time = 02:07:15

Memory = .083 Gb (Passed, expected: 0.2 Gb)

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

Time = 00:00:59

Memory = .084 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 2D_cavityBoussinesq

Run tutorial 2D_cavityBoussinesq ...

Done running tutorial 2D_cavityBoussinesq ...

2D_cavityBoussinesq has converged 


Time = 00:00:46

Memory = .080 Gb (Passed, expected: 0.2 Gb)

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

Time = 00:00:10

Memory = .080 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 2D_FFTF

Run tutorial 2D_FFTF ...

Done running tutorial 2D_FFTF ...

Energy steady state has converged 

Transient has converged 

Reactivity contributions with 0.01 relative error
Perfect match:
|             | Simulated | Expected |
|:------------|:---------:|:--------:|
| Doppler     | 238.1499 | 238.1535 |
| Cladding    | 0.8706511 | 0.8706593 |
| Density     | -2.52138 | -2.523688 |
| Structures  | -5.555345 | -5.55205 |
| Driveline   | -7.253883 | -7.249654 |
| GEM         | -464.8492 | -464.8492 |
| Total power | 28763.92 | 28761.75 |

Time = 02:44:26

Memory = .152 Gb (Passed, expected: 0.3 Gb)

--------------------------------------------------------------------------------

## 2D_fullCoupling

Run tutorial 2D_fullCoupling ...

Done running tutorial 2D_fullCoupling ...

Steady State has converged 

Transient (no driveline) has converged 


Time = 00:00:17

Memory = .082 Gb (Passed, expected: 0.2 Gb)

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

Time = 00:52:05

Memory = .088 Gb (Passed, expected: 0.3 Gb)

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

Time = 00:20:47

Memory = .179 Gb (Passed, expected: 0.3 Gb)

--------------------------------------------------------------------------------

## 2D_onePhaseAndPointKineticsCoupling

Run tutorial 2D_onePhaseAndPointKineticsCoupling ...

Done running tutorial 2D_onePhaseAndPointKineticsCoupling ...

Steady State has converged 

Transient (no driveline) has converged 

Transient (with driveline) has converged 

Transient (with boron) has converged 


Time = 00:00:52

Memory = .082 Gb (Passed, expected: 0.2 Gb)

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

Time = 00:01:33

Memory = .082 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 2D_voidMotionNoPhaseChange

Run tutorial 2D_voidMotionNoPhaseChange ...

Done running tutorial 2D_voidMotionNoPhaseChange ...

2D_voidMotionNoPhaseChange has converged 


Time = 00:01:18

Memory = .083 Gb (Passed, expected: 0.2 Gb)

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

Time = 00:27:24

Memory = .788 Gb (Passed, expected: 1.5 Gb)

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


Time = 02:22:16

Memory = .094 Gb (Passed, expected: 0.2 Gb)

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

Time = 00:11:20

Memory = .175 Gb (Passed, expected: 1.0 Gb)

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

Time = 01:27:03

Memory = 1.076 Gb (Passed, expected: 1.5 Gb)
