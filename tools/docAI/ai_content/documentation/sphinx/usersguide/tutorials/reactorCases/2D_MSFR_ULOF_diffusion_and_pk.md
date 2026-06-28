# 2D MSFR ULOF - Diffusion and Point-Kinetics

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-onePhase-blue.svg)]() [![badge](https://img.shields.io/badge/Neutronics-diffusion-blue.svg)]() [![badge](https://img.shields.io/badge/Multiphysics-picardLoop-orange.svg)]()


## Description

`2D_MSFR_ULOF_diffusion_and_pk` is a 2-D *r-z* model of a Molten Salt Fast Reactor. It solves for neutronics and thermal-hydraulics.

This case is the same as `2D_MSFR`, except that:
- a coupled fluid and energy simulation is run to improve the coupling between energy and velocity
- a loss-of-flow transient is simulated by using both diffusion and point kinetics

The tutorial was prepared by Arnaldo Mattioli (Politecnico di Milano) and revised by C. Fiorina (EPFL).