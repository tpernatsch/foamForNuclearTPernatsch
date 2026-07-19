# Neutronics Diffusion: Slab Reactor with a Reflector

Tags: [![badge](https://img.shields.io/badge/Neutronics-diffusion-blue.svg)]()

## Description

This tutorial is the second in a series of neutronics tutorials. In this tutorial, GeN-Foam is used to demonstrate how a neutron reflector reduces the size of the fissile zone and enables criticality.

<img src="images/2_reactorSlabReflected_1D_1Gr_neutronicDiffusion_results_mesh.png" width="250"/>

*Fig. 1. Slab reactor with reflector visualized with ParaView (core in red, reflector in white).*

## Theory

We start from the diffusion equation in two regions: the Core ($C$) and the Reflector ($R$). The core has length $L$, and the reflector has thickness $d$.

$$
D_C \nabla^2 \phi_C + \nu \Sigma_f^C \phi_C - \Sigma_a^C \phi_C = 0
$$

$$
D_R \nabla^2 \phi_R - \Sigma_a^R \phi_R = 0
$$

The core is centered at $z=0$. Therefore, the solutions to the above equations are

$$
\phi_C(z) = A \cos(B_g z)
$$

$$
\phi_R(z) = E \exp\left(\frac{-z}{L_R}\right) + F \exp\left(\frac{z}{L_R}\right)
$$

where
$L_R = \sqrt{\frac{D_R}{\Sigma_a^R}}$
and
$B_g = \frac{\nu \Sigma_f^C - \Sigma_a^C}{D_C}$.

By imposing zero neutron flux at the outer boundary of the reflector,

$$
\phi_R \left(z=\frac{L}{2}+d \right) = 0
\quad \rightarrow \quad
F = - E \exp \left( -\frac{L+2d}{L_R} \right)
$$

Using continuity of the neutron flux between the core and the reflector, we obtain

$$
\phi_R(z) = A \frac{\cos \left(B_g \frac{L}{2} \right)}{\sinh \left( \frac{d}{L_R} \right)} \sinh \left( \frac{-z+L/2+d}{L_R} \right)
$$

Next, we enforce current continuity at the core–reflector interface to determine $B_g$:

$$
D_C \frac{\partial \phi_C}{\partial z} = D_R \frac{\partial \phi_R}{\partial z}
\quad \rightarrow \quad
D_C B_g \tan \left( B_g \frac{L}{2} \right) = \frac{D_R}{L_R} \coth \left( \frac{d}{L_R} \right)
$$

## Simulation

The case is prepared to verify the analytic results above. You can change the geometry in the [`system/neutroRegion/blockMeshDict`](./system/neutroRegion/blockMeshDict) file using the `lCore` parameter for the core length and `lRefl` for the reflector thickness. The nuclear data are specified in the [`constant/neutroRegion/nuclearData`](./constant/neutroRegion/nuclearData) file.

Before running the case, clean it with the following command:

```bash
./Allclean
```

Then run the case using:

```bash
./Allrun
# or
./Alltest
```

The case should take no more than a few seconds. The results are stored in the `1` folder. The `1/uniform/reactorState` file contains the integral simulation parameters (e.g., $k_\text{eff}$ and the power).

## Analysis

To compare the GeN-Foam results with the analytic solution, run the following command:

```bash
#       Script      timeStep  fieldName
#       v           v         v
python3 Allplot.py  1         flux0
```

![](./images/2_reactorSlabReflected_1D_1Gr_neutronicDiffusion_results_neutronFlux.png)

*Fig. 2. Flux distribution computed by GeN-Foam. Comparison with the analytic solution with and without a reflector.*

From Figure 2, we can see that the reactor requires less fissile material to reach criticality.

## Possible exercises

- How is $k_\text{eff}$ affected when increasing or decreasing the reflector thickness?