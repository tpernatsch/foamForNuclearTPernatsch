# PWR core

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-onePhase-blue.svg)]() [![badge](https://img.shields.io/badge/Neutronics-diffusion-blue.svg)]() [![badge](https://img.shields.io/badge/Multiphysics-looseCoupling-orange.svg)]()


## Description

Simple PWR core model using coupled neutronics diffusion and one-phase thermal-hydraulics sub-solvers. Only one fuel zone is used.


## How to run

```bash
python3 Allclean.py

python3 Allrun.py
```

## Results

<img src="./images/fig_results_fluidRegion_400_T.png" width=500>

*Fig 1. Temperature distribution in 3D.*

<div>
    <img src="./images/fig_results_neutroRegion_400_TCool_yz.png" width=500>
    <img src="./images/fig_results_neutroRegion_400_TClad_yz.png" width=500>
    <img src="./images/fig_results_neutroRegion_400_TFuel_yz.png" width=500>
</div>

*Fig 2. Temperature distributions of coolant, cladding and fuel in XY-projection.*
