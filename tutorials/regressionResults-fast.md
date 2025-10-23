# GeN-Foam regression test suite

Regression test of 29-sept.-2025 - 10:10:49  
OpenFOAM Version: v2506  
Tutorials successfully completed: 11/11  [![regressionTest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]()  

## Summary

| Tutorial | Status | Time | GeN-Foam memory usage [Gb] |
|:---------|:------:|:-:|:-:|
| guidedCases/1_reactorSlab_1D_1Gr_neutronicDiffusion | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:01 | .027 < 0.2 |
| guidedCases/2_reactorSlabReflected_1D_1Gr_neutronicDiffusion | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:01 | .026 < 0.2 |
| featureCases/1D_HX/onePhase | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:06 | .087 < 0.2 |
| featureCases/1D_PSBT_SC/Phase_Ex1_12223 | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:10 | .082 < 0.2 |
| featureCases/2D_cavityBoussinesq | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:31 | .083 < 0.2 |
| featureCases/2D_externalSourceDiffusion | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:02 | .084 < 0.2 |
| featureCases/2D_fullCoupling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:08 | .086 < 0.2 |
| featureCases/2D_onePhaseAndPointKineticsCoupling | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:30 | .085 < 0.2 |
| reactorCases/1D_thermalMSR_pointKinetics | [![Alltest](https://img.shields.io/badge/Alltest-passed-color.svg?style=flat-square)]() | 00:00:43 | .087 < 0.2 |

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

Memory = .027 Gb (Passed, expected: 0.2 Gb)

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

Time = 00:00:01

Memory = .026 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/1D_HX/onePhase

Run tutorial 1D_HX onePhase ...

Done running tutorial 1D_HX onePhase ...

1D_HX onePhase has converged 


Time = 00:00:06

Memory = .087 Gb (Passed, expected: 0.2 Gb)

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

Time = 00:00:10

Memory = .084 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/2D_cavityBoussinesq

Run tutorial 2D_cavityBoussinesq ...

Done running tutorial 2D_cavityBoussinesq ...

2D_cavityBoussinesq has converged 


Time = 00:00:31

Memory = .086 Gb (Passed, expected: 0.2 Gb)

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

Time = 00:00:02

Memory = .086 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/2D_flowOverHeatedPlate

Run tutorial 2D_flowOverHeatedPlate ...

Done running tutorial 2D_flowOverHeatedPlate ...

2D_flowOverHeatedPlate has converged 

Perfect match:
|           | Simulated | Expected |
|:----------|:---------:|:--------:|
| TMax      |  309.143  |  309.143 |

Time = 00:01:20

Memory = .125 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/2D_fullCoupling

Run tutorial 2D_fullCoupling ...

Done running tutorial 2D_fullCoupling ...

Steady State has converged 

Transient (no driveline) has converged 


Time = 00:00:07

Memory = .088 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## featureCases/2D_onePhaseAndPointKineticsCoupling

Run tutorial 2D_onePhaseAndPointKineticsCoupling ...

Done running tutorial 2D_onePhaseAndPointKineticsCoupling ...

Steady State has converged 

Transient (no driveline) has converged 

Transient (with driveline) has converged 

Transient (with boron) has converged 


Time = 00:00:30

Memory = .088 Gb (Passed, expected: 0.2 Gb)

--------------------------------------------------------------------------------

## reactorCases/1D_thermalMSR_pointKinetics

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

Time = 00:00:43

Memory = .089 Gb (Passed, expected: 0.2 Gb)
