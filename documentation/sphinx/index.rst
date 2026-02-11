.. warning::

   This code, tutorials and documentation are currently being finalized and we expect them to be complete by the end of February 2026. This process is expected to lead to some modification in dictionaries and keywords

##############
foamForNuclear
##############

**foamForNuclear** is a **general-purpose, OpenFOAM®-based multiphysics platform**
for the analysis and design of nuclear systems. It was established in 2025 through
the **integration and extension** of two mature OpenFOAM-based projects—**GeN-Foam**
and **OFFBEAT**—both of which had been developed, verified, and validated over
more than a decade. GeN-Foam was originally conceived for reactor-scale
multiphysics simulations, while OFFBEAT focused on advanced nuclear fuel
performance modeling and nonlinear thermo-mechanics.

Building on this combined heritage, **foamForNuclear** provides a
**modular and extensible framework** capable of simulating a broad range of
coupled physics, from core neutronics and thermal-hydraulics to advanced
thermomechanics and detailed fuel behavior. Each physical model is implemented
as an independent **module** that solves a specific set of governing equations.
These modules are accessed through two user-facing applications:

- **GeN-Foam** is the default, general-purpose multiphysics application. It supports loose or tight coupling of an arbitrary number of physics on regions that may be independent, overlapping, or interacting through coupled boundaries.
- **OFFBEAT** is a specialized application dedicated to nuclear fuel behavior, with capabilities covering standard LWR fuel, oxide-based fast-reactor fuel, metallic fuel, and TRISO fuel.

**foamForNuclear** is shipped with a comprehensive Python API for programmatic
pre- and post-processing.


----

**Repository:** |gitlab-badge|

**Doxygen documentation:** |doxygen-badge|

----

.. |gitlab-badge| image:: https://img.shields.io/badge/foamForNuclear-GitLab-orange?logo=gitlab
   :target: https://gitlab.com/foamForNuclear/foamForNuclear
   :alt: foamForNuclear GitLab

.. |doxygen-badge| image:: https://img.shields.io/badge/foamForNuclear-Doxygen-blue?logo=doxygen
   :target: https://foamfornuclear.gitlab.io/foamForNuclear/doxygen/index.html
   :alt: foamForNuclear Doxygen

.. note::

   The master branch contains stable software and is normally
   updated every six months, soon after major OpenFOAM releases. The develop branch contains more recent developments.
   Both branches undergo comprehensive testing before any commit.



.. toctree::
   :hidden:
   :maxdepth: 4

   Home <self>
   releasenotes/index
   usersguide/index
   cppapi/index
   pythonapi/index
   theoryReferences
   V_Vreferences
   contributors
   howToContribute
