# 1D PSBT SC

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-twoPhase-blue.svg)]()

The OECD/NRC PWR PSBT benchmark was organized based on the NUPEC database. It is a well-known benchmark for code validation.

## Phase I Exercise 1 12223

This tutorial presents a simple 1-D case for water boiling based on *OECD/NRC Benchmark based on NUPEC Pressurised Water Reactor (PWR) Subchannel and Bundle Tests (PSBT)*, case 12223 (exercise 1).

## Phase II

This section presents three exercises of *Pressurised Water Reactor (PWR) Subchannel and Bundle Tests (PSBT) - Phase II*.

### Exercise 1 01-5215

Phase II Exercise 1 is in the water boiling condition. The goal is to calculate the liquid temperature at a specific position.

In this exercise, we use similar boundary conditions (defined in `0/fluidRegion`), the same model (defined in `constant/fluidRegion`), and the same solution strategy (defined in `system`).

This exercise uses a challenging power gradient in the radial direction. The model aims to calculate the liquid temperature at an axial height of 4.775 m for comparison with experimental values.

### Exercise 2 04-6770

Phase II Exercise 2 is in the steady-state DNB condition. In these cases, the power level increases gradually toward the critical heat flux condition.

In this exercise, we use similar boundary conditions (defined in `0/fluidRegion`), the same model (defined in `constant/fluidRegion`), and the same solution strategy (defined in `system`).

This exercise uses distributions in both the axial and radial directions. Therefore, we added a `powerDensity.fixedPower` file in `0/fluidRegion` to represent the initial axial power distribution. Then, we defined the power increase ratio relative to the initial value at each time position (see the `powerTimeProfile` block in `constant/fluidRegion/phaseProperties`).

In this exercise, we need to model the transition from single-phase liquid to two-phase water, and then to the critical heat flux condition. Therefore, we implemented a CHF look-up table model to predict the CHF value.

We can obtain the DNB power and the DNB axial position for each case and compare them with experimental values.

### Exercise 3 11-0312

Phase II Exercise 3 is in the transient DNB condition. In these cases, the boundary condition changes over time. Each case also gradually reaches the CHF condition through changes in the boundary condition.

In this exercise, we use different boundary conditions (defined in `0/fluidRegion`), the same model (defined in `constant/fluidRegion`), and the same solution strategy (defined in `system`).

This exercise uses distributions in both the axial and radial directions. Therefore, we added a `powerDensity.fixedPower` file in `0/fluidRegion` to represent the initial axial power distribution. Then, we defined the power increase ratio relative to the initial value at each time position (see the `powerTimeProfile` block in `constant/fluidRegion/phaseProperties`).

This exercise covers transient DNB cases. The boundary condition changes over time.

In this exercise, we also need to model the transition from single-phase liquid to two-phase water, and then to the critical heat flux condition. Therefore, the CHF look-up table is implemented in this exercise as well.

We can obtain the DNB power and the time at which DNB occurs for each case and compare them with experimental values.

## Run simulation

All cases use the same run script. Simply execute the following commands in the terminal.

```bash
./Allclean
./Allrun
# or
./Alltest
```