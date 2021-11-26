# GeN-Foam README file {#README}

GeN-Foam is a multi-physics solver for reactor analysis based on OpenFOAM \cite 2021OpenFOAM \cite Weller1998 (ESI distribution from www.openfoam.com ! currently 2106). It can solve (coupled or alternatively) for:

- neutronics, with models for point kinetics, diffusion (transient and eigenvalue), adjoint diffusion (only eigenvalue), SP3 (transient and eigenvalue), discrete ordinates (only eigenvalue);
- one-phase thermal-hydraulics, according to both RANS-CFD and porous-medium coarse-mesh approaches (the two approaches can be combined in the same mesh);
- two-phase porous-medium thermal-hydraulics, according to an Euler-Euler model, with full capabilities for sodium-cooled fast-reactors and pre-CHF capabilities for ligh water reactors;
- thermal-mechanics based on linear thermo-elasticity, which can be used to evaluate deformations in a core. Such deformations are used to modify the meshes for thermal-hydraulics and neutronics.
A 1-D sub-scale model is also employed cell-by-cell for calculating fuel temperatures in coarse-mesh models of a core.  It should be mentioned that GeN-Foam was mainly designed for coarse-mesh analyses of a reactor core (with porous medium approach and sub-scale representation of fuel), and not for detailed pin-by-pin models.

GeN-Foam is an unusually complex OpenFOAM solver. For this reason, some documentation (associated to the repository, both in the form of a gitLab wiki and of a Doxygen-generated documentation) has been prepared to facilitate its use. In addition, several commented tutorials have been prepared to showcase use and capabilities of the solver. An EMPTY case is also provided that can be used for step-by-step building one’s own case. It is recommended to start from the EMPTY case to build each new case, as it already includes a consistent minimum set of (dummy) files that have to be present independent of the physics that are solved for. Beside this documentation, users ar encouraged to make use of the
typical OpenFOAM resources:
- the above-mentioned user guide and basic tutorials;
- the high-level C++-based object-oriented language of OpenFOAM, which normally allows to easily understand the logic of a solver;
- the support of the community.

Please notice that a new version of OpenFOAM is released twice a year. It may take a few weeks for the developers to update GeN-Foam to a new OpenFOAM release.

N.B.: GeN-Foam is a flexible tool that allows modeling irregular geometries and particularly complex phenomena. However, it requires a good familiarity with OpenFOAM and a very solid back-ground in multi-physics nuclear applications, with particular regard to CFD. In addition, good familiarity with C++ and the OpenFOAM API will be necessary to unlock the full potential of the code. The OpenFOAM API and the class-based structure of GeN-Foam allows an experienced user to quickly an safely add models and equations, thus tailoring the code to their needs.

© All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (EPFL), Switzerland, 2021
- Main author of the code: Carlo Fiorina
- Main author of the thermal-hydraulics class and of the point-kinetics solver: Stefan Radman
- Other contributions are individually acknowledged in the haeader files