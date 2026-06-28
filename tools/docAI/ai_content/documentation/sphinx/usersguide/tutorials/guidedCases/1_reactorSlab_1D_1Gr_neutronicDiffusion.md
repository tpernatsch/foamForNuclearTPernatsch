# Neutronics Diffusion: Slab Reactor

Tags: [![badge](https://img.shields.io/badge/Neutronics-diffusion-blue.svg)]()

## Description

This tutorial is a simple demonstration of GeN-Foam neutronic solvers. We solve the classical slab reactor problem using one energy group and no reflectors.

The one-group neutronics diffusion equation can be written as:

$$
\underbrace{D \nabla^2 \phi}_\text{diffusion} + \underbrace{\nu \Sigma_f \phi}_\text{fission} \underbrace{- \Sigma_a \phi}_\text{absorption} = 0
$$

For a slab reactor of size $L$, the diffusion equation can be simplified to:

$$
D \frac{\partial^2 \phi}{\partial z^2} + \nu \Sigma_f \phi - \Sigma_a \phi = 0
$$

We then introduce the material buckling $B_m^2$ into the diffusion equation:

$$
\frac{\partial^2 \phi}{\partial z^2} + B_m^2 \phi = 0
$$

with:

$$
B_m^2 = \frac{k_\infty - 1}{M^2}
\quad \quad
k_\infty = \frac{\nu \Sigma_f}{\Sigma_a}
\quad \quad
M = \sqrt{\frac{D}{\Sigma_a}}
$$

The solution of this differential equation has the form:

$$
\phi(z) = A \sin \left( B_g z \right)
$$

where $A$ is the normalization constant and $B_g$ is the geometrical buckling.

The boundary conditions imposed to the slab reactor correspond to a flux that tends to zero. Therefore, for the fundamental mode, $B_g$ must be equal to $\frac{\pi}{L}$:

$$
\phi(z) = A \sin \left( \frac{\pi z}{L} \right)
$$

The criticality condition closes the system, which implies:

$$
B_m^2 = B_g^2
\rightarrow
\frac{k_\infty - 1}{M^2} = \left( \frac{\pi}{L} \right)^2
$$

If we assume the following values:
- $\nu \Sigma_f = 5.5675$
- $\Sigma_a = 5.0082$
- $D = 0.1275$

then the resulting $k_\infty$ is $\approx 1.11168$. The reactor should reach criticality for a length $L$ of:

$$
L = \pi \sqrt{\frac{M^2}{k_\infty - 1}}
= \pi \sqrt{\frac{D}{\nu \Sigma_f - \Sigma_a}}
= 1.5 \text{ m}
$$

## Simulation

The case has been prepared to verify the above results. You can change the geometry in the [`system/neutroRegion/blockMeshDict`](./system/neutroRegion/blockMeshDict) file using the `dz` parameter. The nuclear data are set in the [`constant/neutroRegion/nuclearData`](./constant/neutroRegion/nuclearData) file.

Before running the case, clean it using the following command:
```bash
./Allclean
```

Then run the case with:
```bash
./Allrun
# or
./Alltest
```

The case should take no more than a few seconds. The results are stored in the `1` folder. The `1/uniform/reactorState` file stores the integral parameters of the simulation (e.g., $k_\text{eff}$ and the power). The flux distribution is stored in the `1/neutroRegion/flux0` file.

<img src="images/1_reactorSlab_1D_1Gr_neutronicDiffusion_results_mesh_flux0.png" width="250"/>

*Fig 1. Neutron flux distribution in the slab reactor (visualized using [ParaView](https://www.paraview.org/)).*

## Analysis

To compare the GeN-Foam results with the analytic solution, run the following command:

```bash
#       Script      timeStep  fieldName
#       v           v         v
python3 Allplot.py  1         flux0
```

The script should output the flux shape and compare it to the analytic expression.

![](./images/1_reactorSlab_1D_1Gr_neutronicDiffusion_results_neutronFlux.png)

*Fig 2. Neutron flux comparison between GeN-Foam and an analytical solution.*

## Possible exercises

- How is $k_\text{eff}$ affected when increasing or decreasing the fuel length?
- How is $k_\text{eff}$ affected when increasing or decreasing `nuSigmaEff`, `sigmaRemoval`, or `D` in [`constant/neutroRegion/nuclearData`](./constant/neutroRegion/nuclearData)?
- Does the normalization power `pTarget` in [`0.orig/uniform/reactorState`](./0.orig/uniform/reactorState) affect $k_\text{eff}$?
- Should the neutron flux vanish at the top and bottom boundaries?
    - To check this, you can change the neutron flux boundary conditions from `fixedValue` to an albedo boundary condition `albedoSP3`. Uncomment the albedo boundary condition and comment the `fixedValue` boundary condition in [`0.orig/neutroRegion/defaultFlux`](./0.orig/neutroRegion/defaultFlux).
    - Then clean and run the simulation.
    - The flux should no longer tend to zero at the boundaries. You can change the `isAddExtrapolationDistance` flag in the [`Allplot.py`](./Allplot.py) script to visualize the extrapolation limit and obtain the following results.

![](./images/1_reactorSlab_1D_1Gr_neutronicDiffusion_results_neutronFlux_extrapolated.png)

*Fig 3. Neutron flux axial distribution using albedo boundary conditions. Comparison with the analytical solution assuming zero flux at the boundaries.*