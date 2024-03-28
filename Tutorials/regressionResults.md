# GeN-Foam regression test suite

Regression test of 09-Feb-2024 - 21:57:56  
Tutorials successfully completed: 20/23  [![regressionTest](https://img.shields.io/badge/regressionTest-failed-red.svg?style=flat-square)]()  

--------------------------------------------------------------------------------

## 1D_boiling

Run tutorial 1D_boiling ...

Done running tutorial 1D_boiling ...

1D_boiling has converged 


--------------------------------------------------------------------------------

## 1D_CHF/imposedPower

Run tutorial 1D_CHF imposedPower ...

Done running tutorial 1D_CHF imposedPower ...

1D_CHF imposedPower has converged 

Structure temperature at the top of the channel at 90 seconds
Perfect match:
|          | Simulated | Expected |
|:---------|:---------:|:--------:|
| TwallTop | 2618.65 | 2625.18 |

--------------------------------------------------------------------------------

## 1D_CHF/imposedTemperature

Run tutorial 1D_CHF imposedTemperature ...

Done running tutorial 1D_CHF imposedTemperature ...

1D_CHF imposedTemperature has converged 

Heat flux of the structure at the top of the channel at 40 seconds
Perfect match:
|           | Simulated | Expected |
|:----------|:---------:|:--------:|
| Qwall Top | 33427.4 | 33423.1 |

--------------------------------------------------------------------------------

## 1D_HX/onePhase

Run tutorial 1D_HX onePhase ...

Done running tutorial 1D_HX onePhase ...

1D_HX onePhase has converged 


--------------------------------------------------------------------------------

## 1D_HX/twoPhase

Run tutorial 1D_HX twoPhase ...

Done running tutorial 1D_HX twoPhase ...

1D_HX twoPhase has converged 


--------------------------------------------------------------------------------

## 1D_MSR_pointKinetics

Run tutorial 1D_MSR_pointKinetics ...

Done running tutorial 1D_MSR_pointKinetics ...

Steady state has converged 

Transient has converged 

Transient end has converged 

Power at the end of transient
Perfect match:
|       | Simulated | Expected |
|:------|:---------:|:--------:|
| Power | 6.13787e+08  | 6.13794e+08 |

--------------------------------------------------------------------------------

## 1D_PSBT_SC/Phase_Ex1_12223

Run tutorial 1D_PSBT_SC ...

Done running tutorial 1D_PSBT_SC ...

1D_PSBT_SC has converged 

alpha.vapour at the end of the channel
Perfect match:
|              | Simulated | Expected |
|:-------------|:---------:|:--------:|
| alpha vapour | 0.133727 | 0.133727 |

--------------------------------------------------------------------------------

## 1D_PSBT_SC/PhaseII_Ex1_01_5215

Run tutorial PSBT_Rod_Bundle_CHF PhaseII_Ex1_01_5215 ...

Done running tutorial PSBT_Rod_Bundle_CHF PhaseII_Ex1_01_5215 

PSBT_Rod_Bundle_CHF PhaseII_Ex1_01_5215 has converged 

alpha.vapour at the end of the channel
Perfect match:
|                       | Simulated | Expected |
|:----------------------|:---------:|:--------:|
| Max alpha vapour      | 0.123835 | 0.123835 |
| T at max alpha vapour | 620.178 | 620.178 |

--------------------------------------------------------------------------------

## 1D_PSBT_SC/PhaseII_Ex2_04_6770

Run tutorial PSBT_Rod_Bundle_CHF PhaseII_Ex2_04_6770 ...

Done running tutorial PSBT_Rod_Bundle_CHF PhaseII_Ex2_04_6770 

PSBT_Rod_Bundle_CHF PhaseII_Ex2_04_6770 has converged 

alpha.vapour at the end of the channel
Perfect match:
|                       | Simulated | Expected |
|:----------------------|:---------:|:--------:|
| Max alpha vapour      | 0.319074 | 0.319074 |
| T at max alpha vapour | 625.826 | 625.826 |

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

--------------------------------------------------------------------------------

## 2D_cavityBoussinesq

Run tutorial 2D_cavityBoussinesq ...

Done running tutorial 2D_cavityBoussinesq ...

2D_cavityBoussinesq has converged 


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
| Doppler     | 237.049 | 237.0728 |
| Cladding    | 0.8180091 | 0.8198844 |
| Density     | -2.518649 | -2.525069 |
| Structures  | -1.027231 | -1.017622 |
| Driveline   | -8.622036 | -8.631667 |
| GEM         | -464.8496 | -464.8497 |
| Total power | 29493.69 | 29490.49 |

--------------------------------------------------------------------------------

## 2D_fullCoupling

Run tutorial 2D_fullCoupling ...

Done running tutorial 2D_fullCoupling ...

Steady State has converged 

Transient (no driveline) has converged 


--------------------------------------------------------------------------------

## 2D_KNS37-L22

Run tutorial 2D_KNS37-L22 ...

Done running tutorial 2D_KNS37-L22 ...

2D_KNS37-L22 has converged 

Test with 0.2 relative error
Perfect match:
|                   | Simulated | Expected |
|:------------------|:---------:|:--------:|
| Time              | 12.52914919555637 | 12.56424533593758 |
| alpha.vapour      | 0.168632 | 0.1577348 |
| T.activeStructure | 826.679800 | 821.0752 |

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

--------------------------------------------------------------------------------

## 2D_onePhaseAndPointKineticsCoupling

Run tutorial 2D_onePhaseAndPointKineticsCoupling ...

Done running tutorial 2D_onePhaseAndPointKineticsCoupling ...

Steady State has converged 

Transient (no driveline) has converged 

Transient (with driveline) has converged 

Transient (with boron) has converged 


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

--------------------------------------------------------------------------------

## 2D_voidMotionNoPhaseChange

Run tutorial 2D_voidMotionNoPhaseChange ...

Done running tutorial 2D_voidMotionNoPhaseChange ...

2D_voidMotionNoPhaseChange has converged 


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
