What is OFFBEAT
===============

OFFBEAT is an open-source, multi-dimensional (1D, 2D, 3D) nuclear fuel
performance code for the thermo-mechanical analysis of nuclear fuel
and cladding. It covers both reactor operating conditions and
out-of-reactor scenarios.

OFFBEAT is developed as an application built on top of the OpenFOAM® C++
library. More specifically, it relies on a custom extension of OpenFOAM
for solid thermo-mechanics, tailored to the analysis of nuclear materials
and fuel performance phenomena.

As a fuel performance code, OFFBEAT focuses on the continuum-scale
description of fuel and cladding behaviour. Nuclear-specific phenomena
are represented through dedicated material and phenomenological models.
The code is intended for both steady-state and transient analyses,
including operational transients and accident scenarios such as LOCA and
RIA.

Typical applications of OFFBEAT include:

- steady-state base irradiation analyses, including temperature and
  stress evaluation in fuel and cladding,
- analysis of pellet–cladding mechanical interaction, including
  higher-dimensional phenomena such as cladding ridging induced by
  pellet hourglassing,
- investigation of non-axisymmetric effects, such as missing pellet
  surface or pellet eccentricity,
- accidental and transient fuel behaviour studies,
- analysis of hydrogen behaviour in cladding, with or without a liner.

Within these applications, OFFBEAT supports the modelling of a wide range
of fuel behaviour phenomena, including fuel densification, swelling and
relocation, creep, plasticity, and thermal expansion. It also supports
temperature- and burnup-dependent material properties. Fission gas
diffusion in fuel grains and eventual release are modelled through the
embedded open-source 0D SCIANTIX code.

Dedicated boundary conditions are provided to model gap heat transfer and
mechanical contact. OpenFOAM mapping algorithms are used to transfer
fields across the gap, including in the presence of non-conformal
interfaces.

OFFBEAT is developed with an emphasis on modularity and extensibility.
Physical models, material laws, and numerical options are structured to
facilitate the implementation of new physics and the extension of
existing capabilities.