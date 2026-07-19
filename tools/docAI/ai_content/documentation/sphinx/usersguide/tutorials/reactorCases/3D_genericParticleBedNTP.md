# Particle-bed Nuclear Thermal Propulsion Reactor

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-onePhase-blue.svg)]() [![badge](https://img.shields.io/badge/ThermalMechanics-extendedThermoMechanics-blue.svg)]() [![badge](https://img.shields.io/badge/Multiphysics-looseCoupling-orange.svg)]()

Author: Zach Hughes, Texas A&M University (zhughes@tamu.edu). Date: 05.05.2024

## Instructions

- Run `./Allrun` from the top directory.
- This executes a thermal-hydraulic simulation using the Albedo BC-generated power distribution.

## Introduction

This repository contains a transient, coupled thermal-hydraulic/neutronic solution for a particle bed reactor (PBR) core. The power densities are obtained using Serpent. The core was designed for use in a Nuclear Thermal Propulsion (NTP) system. The PBR design is based on the work by Ludewig (1993).

## References

Details of the calculations and results can be found in:

Zachary Hughes, Nolan MacDonald, Cristian Garza, Vishal Patel, Rok Krpan, Carlo Fiorina, Multiphysics modeling and preliminary characterization of particle-bed reactors for nuclear thermal propulsion, *Progress in Nuclear Energy*, Volume 190, 2026, 105953, ISSN 0149-1970, https://doi.org/10.1016/j.pnucene.2025.105953.

Available at: https://www.sciencedirect.com/science/article/pii/S0149197025003518?dgcid=author

Open-access pre-print available at: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5337606