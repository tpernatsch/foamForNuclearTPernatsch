foamForNuclear
==============

foamForNuclear is a community-driven computational framework for simulating
multiphysics phenomena in nuclear systems. It is built on top of the `OpenFOAM
<https://www.openfoam.com/>`_ computational fluid dynamics (CFD) library and
extends its capabilities to model complex nuclear reactor components and safety
scenarios. foamForNuclear supports steady-state and transient simulations,
enabling analysis of fluid flow, heat transfer, and coupled
neutronics/thermal-hydraulic behavior in reactor systems. Models are constructed
using flexible mesh-based geometries and can incorporate detailed material
properties and boundary conditions relevant to nuclear engineering applications.

foamForNuclear leverages OpenFOAM's object-oriented C++ design, providing a
modular structure for implementing nuclear-specific solvers and physics models.
Parallel computing is achieved through OpenFOAM's domain decomposition and
MPI-based parallelism, allowing high-resolution simulations on large-scale HPC
systems.


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

   This project is under active development.


Contents
--------

.. toctree::
   :maxdepth: 1

   Home <self>
   releasenotes/index
   usersguide/index
   pythonapi/index
   cppapi/index
   theoryReferences
   V_Vreferences
   contributors
