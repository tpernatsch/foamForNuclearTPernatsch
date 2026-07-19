# GeN-Foam regression test suite

Regression test of 14-Jul-2026 - 03:59:41  
OpenFOAM Version: v2606  
Tutorials successfully completed: 10/10  [![regressionTest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]()  

## Summary

| Tutorial | Status | Time | GeN-Foam memory usage [Gb] |
|:---------|:------:|:-:|:-:|
| ../tests | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:01 | 0 < 0.2 |
| guidedCases/1_reactorSlab_1D_1Gr_neutronicDiffusion | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:01 | 0 < 0.2 |
| guidedCases/2_reactorSlabReflected_1D_1Gr_neutronicDiffusion | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:01 | 0 < 0.2 |
| featureCases/1D_HX/onePhase | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:14 | .090 < 0.2 |
| featureCases/1D_PSBT_SC/Phase_Ex1_12223 | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:20 | .085 < 0.2 |
| featureCases/2D_cavityBoussinesq | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:58 | .087 < 0.2 |
| featureCases/2D_externalSourceDiffusion | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:06 | .086 < 0.2 |
| featureCases/2D_fullCoupling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:17 | .089 < 0.2 |
| featureCases/2D_onePhaseAndPointKineticsCoupling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:58 | .089 < 0.2 |
| reactorCases/1D_thermalMSR_pointKinetics | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:02:36 | .090 < 0.2 |

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

All tests passed

Time = 00:00:01

Memory = 0 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## guidedCases/2_reactorSlabReflected_1D_1Gr_neutronicDiffusion

Run tutorial 2_reactorSlabReflected_1D_1Gr_neutronicDiffusion ...

Done running tutorial 2_reactorSlabReflected_1D_1Gr_neutronicDiffusion ...

2_reactorSlabReflected_1D_1Gr_neutronicDiffusion has converged 

All tests passed

Time = 00:00:01

Memory = 0 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/1D_HX/onePhase

Run tutorial 1D_HX onePhase ...

Done running tutorial 1D_HX onePhase ...

1D_HX onePhase has converged 


Time = 00:00:14

Memory = .090 Gb (Passed, expected: 0.2 Gb)

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

Time = 00:00:20

Memory = .085 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/2D_cavityBoussinesq

Run tutorial 2D_cavityBoussinesq ...

Done running tutorial 2D_cavityBoussinesq ...

2D_cavityBoussinesq has converged 


Time = 00:00:58

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

Time = 00:00:06

Memory = .086 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/2D_fullCoupling

Run tutorial 2D_fullCoupling ...

Done running tutorial 2D_fullCoupling ...

Steady State has converged 

Transient (no driveline) has converged 


Time = 00:00:17

Memory = .089 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/2D_onePhaseAndPointKineticsCoupling

Run tutorial 2D_onePhaseAndPointKineticsCoupling ...

Done running tutorial 2D_onePhaseAndPointKineticsCoupling ...

Steady State has converged 

Transient (no driveline) has converged 

Transient (with driveline) has converged 

Transient (with boron) has converged 


Time = 00:00:58

Memory = .089 Gb (Passed, expected: 0.2 Gb)

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

Time = 00:02:36

Memory = .090 Gb (Passed, expected: 0.2 Gb)
