General modelling framework
===========================

OFFBEAT is developed using standard object-oriented design principles and a modular
architecture. Its overall structure is designed to separate physical responsibilities,
numerical solution strategies, and phenomenological modelling. At the same time, it
allows these components to interact within a unified simulation workflow.

The framework is intentionally generic and extensible. This enables the introduction of new
physics, material models, or numerical formulations without requiring (in most cases)
intrusive changes to the core solver structure.


High-level solver organisation
------------------------------

At its core, OFFBEAT addresses the thermo-mechanical behaviour of nuclear fuel and cladding
by solving the heat conduction and momentum balance equations on a shared computational
domain. These governing equations are discretized using the finite volume method on
unstructured meshes of arbitrary geometry.

The thermo-mechanical problem is decomposed into distinct solver components. Each component
is responsible for a specific physical domain (e.g. thermal, mechanical, gap, and plenum
behaviour). These components operate on the same mesh and exchange information through
well-defined interfaces.

The overall solution strategy follows a segregated approach. In this approach, individual
physical problems are solved sequentially within a time step, and their interactions are
resolved through iterative coupling at the top level.


Separation of physics, materials, and phenomenology
---------------------------------------------------

A key aspect of the OFFBEAT framework is the clear separation between:

- governing equations and numerical solvers,
- material constitutive laws,
- irradiation-driven and fuel-specific phenomenological models,
- boundary and interface conditions.

Material behaviour is described through dedicated constitutive models, accounting for
elastic, plastic, and time-dependent effects as appropriate. Nuclear-specific phenomena
such as swelling, relocation, fission gas release, and gap behaviour are handled through
specialised models. These models interact with the thermal and mechanical solvers without
being hard-wired into them.

This separation allows individual models to be developed, verified, replaced, or extended
independently. It also supports different levels of modelling detail depending on the
objectives of a given analysis.


Boundary conditions and interfaces
----------------------------------

OFFBEAT leverages standard OpenFOAM boundary condition mechanisms where applicable. For
example, these mechanisms can be used to prescribe field values or fluxes at domain
boundaries. In addition, dedicated boundary and interface models have been developed to
address fuel performance-specific requirements. These include heat transfer across the
fuel–cladding gap and mechanical contact between non-conformal surfaces.

Mapping techniques provided by the OpenFOAM framework are used to transfer fields across
interfaces. This allows fuel and cladding regions to be discretized independently while
remaining thermo-mechanically coupled.


Development and verification philosophy
---------------------------------------

The development of OFFBEAT proceeds in parallel with verification and validation activities.
Verification cases, separate-effect tests, and comparisons against experimental databases are
used to assess the behaviour of individual models and solver components. This
verification-oriented approach is intended to support both method development and the
progressive extension of the code’s validated application domain.