# **GeN-Foam Documentation**

GeN-Foam is a multi-physics solver for reactor analysis [@FIORINA201524]. It can solve (coupled or alternatively) for:

- neutronics, with models for point kinetics, diffusion (transient and eigenvalue), adjoint diffusion (only eigenvalue), SP3 (transient and eigenvalue), discrete ordinates (only eigenvalue);
- one-phase thermal-hydraulics, according to both RANS-CFD and porous-medium coarse-mesh approaches (the two approaches can be combined in the same mesh);
- two-phase porous-medium thermal-hydraulics, according to an Euler-Euler model, with models and correlations available for sodium and water;
- temperatures in sub-scale solid structures in porous-medium regions, based on user-selectable models including 1-D fuel, fixed temperature, fixed power, fuel pebbles, heated rods, as well as a generic lumped-parameters model based on the concept of electric equivalence;
- thermal-mechanics based on linear thermo-elasticity, which can be used to evaluate deformations and temperatures in solid structures. Deformations can be used to modify the meshes for thermal-hydraulics and neutronics.

It should be mentioned that GeN-Foam was mainly designed for coarse-mesh analyses of a reactor core (with porous medium approach and sub-scale representation of fuel). However, pin-by-pin and other heterogeneous models can be obtained by connecting the thermal-mechanics and thermal-hydraulics regions using coupled boundary conditions.


## Content

- [Installation](installation.md)
- [Preprocessing](user_guide/pre_processing.md)
- [Running GeN-Foam](user_guide/running_GeN-Foam.md)
- [Postprocessing](user_guide/post_processing.md)
- [Tips and tricks](user_guide/tips_and_tricks.md)
- [Important notes](user_guide/important_notes.md)

- [User manual](user_guide/index.md)
	- [Neutronics](user_guide/neutronics.md)
	- [Thermal-hydraulics](user_guide/thermalHydraulics.md)
	- [Thermal-mechanics](user_guide/thermoMechanics.md)
	- [Coupling and time stepping](user_guide/coupling.md)
	- [FMU](user_guide/fmu.md)

- [Tutorials](tutorials/index.md)

- [Recent changes in the case folder](release_notes/index.md)

- [Online Doxygen-generated documentation](https://foam-for-nuclear.gitlab.io/GeN-Foam/doxygen/index.html)


N.B.: While the documentation reported on this page is preliminary but mostly finished, several classes have been created in the last 8 years that still lack proper documentation in the header files. We are working on it...
