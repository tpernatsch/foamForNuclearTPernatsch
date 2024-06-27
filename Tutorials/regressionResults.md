# GeN-Foam regression test suite

Regression test of 25-juin-2024 - 14:02:52  
Tutorials successfully completed: 0/23  [![regressionTest](https://img.shields.io/badge/Alltest-pending-grey.svg?style=flat-square)]()  

## Summary

| Tutorial | Status | Time | Memory usage [Gb] |
|:---------|:------:|:-:|:-:|
| 1D_boiling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 2:12.28 | .084 < 0.2 |
| 1D_CHF/imposedPower | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 28:17.85 | .083 < 0.2 |
| 1D_CHF/imposedTemperature | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 3:44.50 | .084 < 0.2 |
| 1D_HX/onePhase | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 0:06.79 | .091 < 0.2 |
| 1D_HX/twoPhase | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 1:10.26 | .092 < 0.2 |
| 1D_MSR_pointKinetics | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 1:32.69 | .091 < 0.2 |
| 1D_PSBT_SC/Phase_Ex1_12223 | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 0:11.07 | .085 < 0.2 |
| 1D_PSBT_SC/PhaseII_Ex1_01_5215 | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 1:28.53 | .092 < 0.2 |
| 1D_PSBT_SC/PhaseII_Ex2_04_6770 | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 2:01:34 | .088 < 0.2 |
| 1D_thermalMSR_pointKinetics | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 0:52.33 | .090 < 0.2 |
| 2D_cavityBoussinesq | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 0:35.74 | .085 < 0.2 |
| 2D_externalSourceDiffusion | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 0:12.13 | .090 < 0.2 |
| 2D_FFTF | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 4:12:12 | .280 < 0.3 |
| 2D_fullCoupling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 0:16.98 | .085 < 0.2 |
| 2D_KNS37-L22 | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 45:37.74 | .095 < 0.2 |
| 2D_MSFR | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 19:11.41 | .200 < 0.3 |
| 2D_onePhaseAndPointKineticsCoupling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 0:58.49 | .088 < 0.2 |
| 2D_onePhaseAndSubcriticalPointKineticsCoupling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 1:40.23 | .087 < 0.2 |
| 2D_voidMotionNoPhaseChange | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 0:59.72 | .089 < 0.2 |
| 3D_gFHR | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 26:37.49 | .895 < 1.5 |
| 3D_NTPfuelAssembly | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 2:20:48 | .101 < 0.2 |
| 3D_SmallESFR | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 3:58.98 | .139 < 0.3 |
| Godiva_SN | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 1:30:17 | 1.061 < 1.5 |

--------------------------------------------------------------------------------

## 1D_boiling

Run tutorial 1D_boiling ...

Done running tutorial 1D_boiling ...

1D_boiling has converged 


Time = 2:12.28

Memory = .084 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 1D_CHF/imposedPower

Run tutorial 1D_CHF imposedPower ...

Done running tutorial 1D_CHF imposedPower ...

1D_CHF imposedPower has converged 

Structure temperature at the top of the channel at 90 seconds with 0.01 relative error
Perfect match:
|          | Simulated | Expected |
|:---------|:---------:|:--------:|
| TwallTop | 2617.31 | 2617.31 |

Time = 28:17.85

Memory = .083 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 1D_CHF/imposedTemperature

Run tutorial 1D_CHF imposedTemperature ...

Done running tutorial 1D_CHF imposedTemperature ...

1D_CHF imposedTemperature has converged 

Heat flux of the structure at the top of the channel at 40 seconds with 0.01 relative error
Perfect match:
|           | Simulated | Expected |
|:----------|:---------:|:--------:|
| Qwall Top | 33429.8 | 33429.8 |

Time = 3:44.50

Memory = .084 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 1D_HX/onePhase

Run tutorial 1D_HX onePhase ...

Done running tutorial 1D_HX onePhase ...

1D_HX onePhase has converged 


Time = 0:06.79

Memory = .091 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 1D_HX/twoPhase

Run tutorial 1D_HX twoPhase ...

Done running tutorial 1D_HX twoPhase ...

1D_HX twoPhase has converged 


Time = 1:10.26

Memory = .092 Gb (Passed, expected: 0.2 Gb)

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
| Power | 613787000 | 613794000 |

Time = 1:32.69

Memory = .091 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 1D_PSBT_SC/Phase_Ex1_12223

Run tutorial 1D_PSBT_SC ...

Done running tutorial 1D_PSBT_SC ...

1D_PSBT_SC has converged 

alpha.vapour at the end of the channel with 0.001 relative error
Perfect match:
|              | Simulated | Expected |
|:-------------|:---------:|:--------:|
| alpha vapour | 0.131886 | 0.131886 |

Time = 0:11.07

Memory = .085 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 1D_PSBT_SC/PhaseII_Ex1_01_5215

Run tutorial PSBT_Rod_Bundle_CHF PhaseII_Ex1_01_5215 ...

Done running tutorial PSBT_Rod_Bundle_CHF PhaseII_Ex1_01_5215 

PSBT_Rod_Bundle_CHF PhaseII_Ex1_01_5215 has converged 

alpha.vapour at the end of the channel with 0.01 relative error
Perfect match:
|                       | Simulated | Expected |
|:----------------------|:---------:|:--------:|
| Max alpha vapour      | 0.123641 | 0.123835 |
| T at max alpha vapour | 620.18 | 620.178 |

Time = 1:28.53

Memory = .092 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 1D_PSBT_SC/PhaseII_Ex2_04_6770

Run tutorial PSBT_Rod_Bundle_CHF PhaseII_Ex2_04_6770 ...

Done running tutorial PSBT_Rod_Bundle_CHF PhaseII_Ex2_04_6770 

PSBT_Rod_Bundle_CHF PhaseII_Ex2_04_6770 has converged 

alpha.vapour at the end of the channel with 0.01 relative error
Perfect match:
|                       | Simulated | Expected |
|:----------------------|:---------:|:--------:|
| Max alpha vapour      | 0.318288 | 0.319074 |
| T at max alpha vapour | 625.825 | 625.826 |

Time = 2:01:34

Memory = .088 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 1D_thermalMSR_pointKinetics

Run tutorial 1D_MSR_pointKinetics ...

Done running tutorial 1D_MSR_pointKinetics ...

Steady state has converged 

Transient has converged 

Transient end has converged 

Power
Perfect match:
|       | Simulated | Expected |
|:------|:---------:|:--------:|
| Power | 4.15258e+08  | 4.15258e+08 |

Time = 0:52.33

Memory = .090 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 2D_cavityBoussinesq

Run tutorial 2D_cavityBoussinesq ...

Done running tutorial 2D_cavityBoussinesq ...

2D_cavityBoussinesq has converged 


Time = 0:35.74

Memory = .085 Gb (Passed, expected: 0.2 Gb)

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

Time = 0:12.13

Memory = .090 Gb (Passed, expected: 0.2 Gb)

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
| Doppler     | 237.0728 | 237.0728 |
| Cladding    | 0.8198844 | 0.8198844 |
| Density     | -2.525069 | -2.525069 |
| Structures  | -1.017622 | -1.017622 |
| Driveline   | -8.631667 | -8.631667 |
| GEM         | -464.8497 | -464.8497 |
| Total power | 29490.49 | 29490.49 |

Time = 4:12:12

Memory = .280 Gb (Passed, expected: 0.3 Gb)

--------------------------------------------------------------------------------

## 2D_fullCoupling

Run tutorial 2D_fullCoupling ...

Done running tutorial 2D_fullCoupling ...

Steady State has converged 

Transient (no driveline) has converged 


Time = 0:16.98

Memory = .085 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 2D_KNS37-L22

Run tutorial 2D_KNS37-L22 ...

Done running tutorial 2D_KNS37-L22 ...

2D_KNS37-L22 has converged 

Test with 0.2 relative error
Perfect match:
|                   | Simulated | Expected |
|:------------------|:---------:|:--------:|
| Time              | 12.447433891622175 | 12.514 |
| alpha.vapour      | 0.167371 | 0.178709 |
| T.activeStructure | 828.185800 | 830.2673 |

Time = 45:37.74

Memory = .095 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 2D_MSFR

Run tutorial 2D_MSFR ...

Done running tutorial 2D_MSFR ...

Thermal-hydraulic steady state has converged 

Energy steady state has converged 

Transient has converged 

Test keff and power  
Perfect match:
|             | Simulated | Expected |
|:------------|:---------:|:--------:|
| keff        | 0.960287 | 0.960283 |
| Total power | 20000000 | 20000000 |

Time = 19:11.41

Memory = .200 Gb (Passed, expected: 0.3 Gb)

--------------------------------------------------------------------------------

## 2D_onePhaseAndPointKineticsCoupling

Run tutorial 2D_onePhaseAndPointKineticsCoupling ...

Done running tutorial 2D_onePhaseAndPointKineticsCoupling ...

Steady State has converged 

Transient (no driveline) has converged 

Transient (with driveline) has converged 

Transient (with boron) has converged 


Time = 0:58.49

Memory = .088 Gb (Passed, expected: 0.2 Gb)

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

Time = 1:40.23

Memory = .087 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 2D_voidMotionNoPhaseChange

Run tutorial 2D_voidMotionNoPhaseChange ...

Done running tutorial 2D_voidMotionNoPhaseChange ...

2D_voidMotionNoPhaseChange has converged 


Time = 0:59.72

Memory = .089 Gb (Passed, expected: 0.2 Gb)

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

Time = 26:37.49

Memory = .895 Gb (Passed, expected: 1.5 Gb)

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
| keff                    | 1.3194369 | 1.3194342 |   -.00020 | # Passed
| Total power neutro [MW] |   937.000 |       937 |         0 | # Passed
| Total power fluid  [MW] |   937.001 |       937 |   -.00010 | # Passed
| Total mass flow  [kg/s] | 30.54932000 |      31.0 |      1.45 | # Passed
| Outlet temperature  [K] |  1929.981 |      1972 |      2.13 | # Passed

All test passed


Time = 2:20:48

Memory = .101 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## 3D_SmallESFR

Run tutorial 3D_SmallMSFR ...

Done running tutorial 3D_SmallMSFR ...

Steady state has converged 

Transient has converged 

keff and power test with 0.001 relative error
Perfect match:
|             | Simulated | Expected |
|:------------|:---------:|:--------:|
| keff        | 0.936873 | 0.936827 |
| Total power | 799894700 | 800000000 |

Time = 3:58.98

Memory = .139 Gb (Passed, expected: 0.3 Gb)

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

Time = 1:30:17

Memory = 1.061 Gb (Passed, expected: 1.5 Gb)
