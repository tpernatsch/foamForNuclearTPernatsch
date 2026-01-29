Key characteristics of OFFBEAT
==============================

This section describes the main methodological and architectural characteristics that
differentiate OFFBEAT from traditional nuclear fuel performance codes.


OpenFOAM as development framework
---------------------------------

OFFBEAT is developed within the OpenFOAM® C++ framework. While OpenFOAM is widely known for its
use in computational fluid dynamics, it is fundamentally a general-purpose toolbox for the
solution of systems of partial differential equations and for the development of multi-physics
simulation tools.

The choice of OpenFOAM as a development platform provides a number of features that are
particularly relevant for fuel performance analysis, especially for multi-dimensional and
high-fidelity studies, including:

- support for unstructured meshes and complex geometries,
- native handling of multiple regions and fields,
- mapping capabilities for non-conformal interfaces,
- scalable linear solvers and parallel execution,
- a flexible architecture for extending and combining physical models.

These features make OpenFOAM well suited for addressing the coupled thermo-mechanical problem
posed by nuclear fuel behaviour, while avoiding assumptions on geometry, dimensionality, or
problem complexity at the framework level.


Finite volume method for solid thermo-mechanics
-----------------------------------------------

A distinctive characteristic of OFFBEAT is the use of the finite volume method (FVM) for the
solution of solid thermo-mechanical problems. While FVM is widely adopted in computational fluid
dynamics, most nuclear fuel performance codes rely instead on the finite element method (FEM)
or, in some cases, on finite difference approaches.

From a conceptual point of view, FVM and FEM share several similarities and can lead to
equivalent discretizations in simple cases. There is therefore no fundamental separation
between the two approaches in terms of applicability to continuum mechanics problems.
Comparative studies have shown that neither method exhibits a systematic advantage in terms of
accuracy or performance for solid mechanics applications :cite:`Cardiff30`.

In OFFBEAT, the use of FVM offers a number of practical advantages:

- inherent local and global conservation properties,
- memory-efficient and scalable iterative solution strategies,
- straightforward formulation based on control-volume balances,
- natural compatibility with other finite-volume-based physics.

These aspects are particularly advantageous in the context of large-scale, multi-dimensional,
and multi-physics fuel performance simulations. The relative simplicity of the finite-volume
formulation is also an important factor, especially in a collaborative development environment
involving contributors with diverse backgrounds in physics and engineering rather than in
numerical method development.


Native multi-dimensional capability
-----------------------------------

OFFBEAT is conceived as a fully multi-dimensional fuel performance code. The underlying
OpenFOAM framework is inherently three-dimensional, and lower-dimensional analyses are
obtained as special cases through appropriate geometrical representations and boundary
conditions.

As a consequence, OFFBEAT is not tied to a predefined geometry, dimensional reduction, or
canonical discretization. Fuel behaviour problems can be formulated using one-, two-, or
three-dimensional models depending on the level of detail required, while relying on the same
solver infrastructure and physical models.

The use of unstructured finite-volume meshes allows complex and non-axisymmetric geometries to
be addressed without imposing restrictions on mesh topology or coordinate system. This makes
it possible to investigate local and asymmetric effects that cannot be captured within
traditional reduced-dimensional fuel performance approaches.

The same modelling workflow applies across spatial dimensionalities, enabling a consistent
transition from engineering-scale analyses to higher-fidelity multi-dimensional simulations
without changes in the underlying numerical framework.


Multi-physics and coupling philosophy
-------------------------------------

Fuel behaviour problems are intrinsically multi-physics, involving the interaction of thermal,
mechanical, and irradiation-driven phenomena. OFFBEAT addresses this complexity by adopting
a modular and segregated solution philosophy. Individual physical problems are formulated and
solved as separate components, while their interactions are handled through iterative coupling
at the solver level.

This design prioritizes clarity, extensibility, and robustness. Physical models can be
introduced, modified, or replaced independently, without requiring a complete reformulation of
the numerical system.

Rather than enforcing a fixed coupling strategy or a predefined level of fidelity, OFFBEAT
supports different modelling choices depending on the problem being addressed. Individual
phenomena can be activated or deactivated according to the objectives of a given analysis,
allowing OFFBEAT to focus on specific mechanisms of interest without introducing unnecessary
model complexity.
