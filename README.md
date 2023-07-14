# GeN-Foam README file {#README}

GeN-Foam is a multi-physics solver for reactor analysis based on OpenFOAM (ESI/OpenCFD distribution from [www.openfoam.com](https://www.openfoam.com), currently v2306). It can solve (coupled or alternatively) for:

- neutronics, with models for point kinetics, diffusion (transient and eigenvalue), adjoint diffusion (only eigenvalue), SP3 (transient and eigenvalue), discrete ordinates (only eigenvalue);
- one-phase thermal-hydraulics, according to both RANS-CFD and porous-medium coarse-mesh approaches (the two approaches can be combined in the same mesh);
- two-phase porous-medium thermal-hydraulics, according to an Euler-Euler model, with full capabilities for sodium-cooled fast-reactors and pre-CHF capabilities for light water reactors;
- thermal-mechanics based on linear thermo-elasticity, which can be used to evaluate deformations in a core. Such deformations are used to modify the meshes for thermal-hydraulics and neutronics.
A 1-D sub-scale model is also employed cell-by-cell for calculating fuel temperatures in coarse-mesh models of a core.  It should be mentioned that GeN-Foam was mainly designed for coarse-mesh analyses of a reactor core (with porous medium approach and sub-scale representation of fuel), and not for detailed pin-by-pin models.

GeN-Foam is an unusually complex OpenFOAM solver. For this reason, some documentation (in the form of an [online Doxygen-generated documentation](https://foam-for-nuclear.gitlab.io/GeN-Foam/index.html) has been prepared to facilitate its use. In addition, the slides of introductory lectures to both OpenFOAM and GeN-Foam are provided in the folder *Documentation/someUsefulDocumentsAndPResentations*. These lectures are taken from an IAEA e-learning course available at https://elearning.iaea.org/m2/course/view.php?id=1286. The course requires registration and a NUCLEUS account, but it should be available to all IAEA member states. Finally, several commented tutorials have been prepared to showcase the use and capabilities of the solver. An EMPTY case is also provided that can be used for step-by-step building one’s own case. It is recommended to start from the EMPTY case to build each new case, as it already includes a consistent minimum set of (dummy) files that must be present independent of the physics that are solved for. Besides this documentation, users are encouraged to make use of the typical OpenFOAM ways:
- the high-level C++-based object-oriented language of OpenFOAM, which normally allows understanding the logic of a solver easily;
- the comments that are typically available in the source code and, in particular, in the header files of each class;
- the support of the community.

Please notice that a new version of OpenFOAM is released by ESI/OpenCFD twice a year. It may take a few weeks for the developers to update GeN-Foam to a new OpenFOAM release.

N.B.: GeN-Foam is a flexible tool that allows the modeling of irregular geometries and particularly complex phenomena. However, it requires a good familiarity with Linux and OpenFOAM, as well as a solid background in multi-physics nuclear applications. Familiarity with CFD methods is strongly recommended. In addition, a good familiarity with C++ and the OpenFOAM API will be important to unlock the full potential of the code. The OpenFOAM API and the class-based structure of GeN-Foam allow an experienced user to quickly and safely add solvers, models, and equations, thus tailoring the code to their needs.

© Contributions are individually acknowledged in the header files
