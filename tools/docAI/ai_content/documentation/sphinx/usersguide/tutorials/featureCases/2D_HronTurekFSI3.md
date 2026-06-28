# Hron–Turek FSI3 Benchmark

Tags: [![badge](https://img.shields.io/badge/ThermalHydraulics-rhoPimpleFoam-blue.svg)]() [![badge](https://img.shields.io/badge/ThermalMechanics-extendedThermoMechanics-blue.svg)]() [![badge](https://img.shields.io/badge/Multiphysics-fsiLoop-orange.svg)]()


## Overview

The Hron–Turek benchmark is a classical 2-D fluid–structure interaction (FSI) test. It consists of laminar incompressible flow in a channel past a rigid circular cylinder, with an elastic rectangular beam attached to the downstream side of the cylinder.  

The **FSI3** variant is the most demanding of the series. It targets strong added-mass coupling and large-amplitude, self-excited oscillations of the flexible beam. Therefore, it is a standard test to validate both partitioned and monolithic FSI solvers.


## Geometry

- **Channel length:** $L = 2.5\ \text{m}$
- **Channel height:** $H = 0.41\ \text{m}$
- **Cylinder center (from left/bottom corner):** $(0.2,\ 0.2)\ \text{m}$
- **Cylinder radius:** $r = 0.05\ \text{m}$
- **Elastic beam (attached to cylinder):**
  - Length: $l = 0.35\ \text{m}$
  - Thickness: $h = 0.02\ \text{m}$
  - Right-bottom corner of undeformed beam: $(0.6,\ 0.19)\ \text{m}$
- **Reference point on beam tip:** $A(0) = (0.6,\ 0.2)\ \text{m}$

The geometry was generated using `blockMesh`. The `blockMeshDict` files were taken from the public repository of solids4Foam: [https://www.solids4foam.com/tutorials/more-tutorials/fluid-solid-interaction/HronTurekFsi3.html](https://www.solids4foam.com/tutorials/more-tutorials/fluid-solid-interaction/HronTurekFsi3.html)

![FSI3 Geometry](images/2D_HronTurekFSI3_problemGeometry.png)


## Governing Equations

- **Fluid:** incompressible Newtonian (Navier–Stokes, laminar)
- **Solid:** compressible elastic, Saint-Venant–Kirchhoff model
- **FSI interface:** no-slip (velocity continuity) and traction equilibrium


## Boundary & Initial Conditions

- **Inlet ($x = 0$):** parabolic profile with mean velocity $\bar{U}$ and maximum $1.5\bar{U}$:
  $$
  v_{\text{in}}(y) = 1.5\,\bar{U}\ \frac{y(H-y)}{(H/2)^2}
  $$
  Ramp smoothly from zero during the first 2 seconds:
  $$
  v(t) = v_{\text{in}} \cdot \frac{1 - \cos\left(\frac{\pi}{2}t\right)}{2}, \quad t < 2\ \text{s}
  $$
- **Walls (top/bottom):** no-slip
- **Cylinder surface:** rigid, no-slip
- **Beam root:** clamped to cylinder
- **Outlet ($x = L$):** zero-mean pressure (reference pressure at outlet)


## Material Properties

|            Parameter             | Value (FSI3) |   Units    |
| :------------------------------: | :----------: | :--------: |
|    Fluid density, $\rho_f$     |     1000     | kg/m$^3$ |
|    Fluid viscosity, $\nu_f$    |    0.001     | m$^2$/2  |
| Mean inlet velocity, $\bar{u}$ |      2       |    m/s     |
|    Solid density, $\rho_s$     |     1000     | kg/m$^3$ |
|  Solid Young's modulus, $E_s$  |     5.6      |    MPa     |
| Solid Poisson's ratio, $\nu_s$ |     0.4      |            |


## Numerical Guidelines

- **Inlet ramp:** use a cosine ramp for the first 2 s to avoid transients
- **Time step:** $\Delta t \approx 10^{-3}\ \text{s}$ (refine if needed)
- **Mesh:** refine near the cylinder/beam interface
- **Coupling:** for partitioned solvers, use acceleration/relaxation to handle the strong added-mass effect


## Expected Results

The benchmark parameters used to validate the results are the beam tip horizontal and vertical displacement histories, as well as the dominant oscillation frequency.

|           Parameter           |    foamForNuclear      | Expected Results [1] |
| :---------------------------: | :--------------------: | :------------------: |
| $u_x$ in mm |  −2.67 ± 2.59 [10.81]  | −2.69 ± 2.53 [10.9]  |
| $u_y$ in mm |   1.53 ± 33.45 [5.38]   |  1.48 ± 34.38 [5.3]  |


## How to run

To run the case, launch the `Allrun` script with `./Allrun`. Note that, due to the large deformations of the solid body, convergence is relatively slow and requires strong under-relaxation of the FSI coupling (0.05). The case takes about 15 hours to converge. Domain decomposition does not help because the mesh size is relatively small.  

The solid body oscillations cease to vary after 4 seconds. Therefore, the computation time can be reduced by simulating until 4 seconds instead of 6 seconds. Future work to accelerate the convergence of tightly coupled FSI problems like this one is foreseen.


## Reference

S. Turek & J. Hron — *Proposal for numerical benchmarking of fluid–structure interaction between an elastic object and laminar incompressible flow* (2006)