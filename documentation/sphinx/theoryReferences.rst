.. _theory_references:

=================
Theory references
=================

--------------
Code structure
--------------

.. _FIORINA201524:

- C. Fiorina, I. Clifford, M. Aufiero, and K. Mikityuk, "`GeN-Foam: A novel
  OpenFOAM® based multi-physics solver for 2D/3D transient analysis of nuclear
  reactors <https://doi.org/10.1016/j.nucengdes.2015.05.035>`_",
  Nuclear Engineering and Design, vol. 294, pp. 24-37, Dec. 2015.

  - Relevance - **GeN-Foam**. The first journal publication on GeN-Foam. It provides a
    general overview of the solver and represents a solid basis to understand the software
    functioning and philosophy.
  - Warning - Much of the structure of the solver has changed.
    The new structure is well described in :ref:`Nervi et al. (2026) <NERVI2026106250>`. Notable changes are:
    generalized coupling loops and generalized coupling with multiple physics and multiple
    geometries (overlapped or not).


.. _SCOLARO2020110416:

- A. Scolaro, I. Clifford, C. Fiorina, and A. Pautz, "`The OFFBEAT multi-dimensional
  fuel behavior solver <https://doi.org/10.1016/j.nucengdes.2019.110416>`_",
  Nuclear Engineering and Design, vol. 358, 110416 2020.

  - Relevance - **OFFBEAT**. The first journal publication on OFFBEAT. @Ale
  - Warning - @Ale


.. _NERVI2026106250:

- G. Nervi, T. Guilbaud, C. Fiorina, M. Hursin, and A. Scolaro, "`foamForNuclear:
  A unified OpenFOAM multi-physics platform for nuclear applications
  <https://doi.org/10.1016/j.pnucene.2026.106250>`_", Progress in Nuclear Energy,
  vol. 193, 2026, 106250, ISSN 0149-1970.


-------------------------------------------
Thermal-hydraulics / heat and mass transfer
-------------------------------------------

- C. Fiorina, I. Clifford, M. Aufiero, and K. Mikityuk, "`GeN-Foam: A novel
  OpenFOAM® based multi-physics solver for 2D/3D transient analysis of nuclear
  reactors <https://doi.org/10.1016/j.nucengdes.2015.05.035>`_",
  Nuclear Engineering and Design, vol. 294, pp. 24-37, Dec. 2015.

  - Relevance - **GeN-Foam**. Good reference for hybrid porous and clear fluid CFD,
    sub-scale sources
    and sinks,  the concept of sub-scale structures, and turbulence in porous media.
    It briefly explains the concept of global and local coordinate systems.
  - Warning - Navier-Stokes equations are now solved based on the real fluid velocity,
    not Darcy velocity (See :ref:`RADMAN2021111178 <RADMAN2021111178>`).

.. _RADMAN2021111178:

- S. Radman, C. Fiorina, and A. Pautz, "`Development of a novel two-phase flow solver
  for nuclear reactor analysis: algorithms, verification and implementation in OpenFOAM
  <https://doi.org/10.1016/j.nucengdes.2021.111178`_", Nuclear Engineering and
  Design, vol. 379, 111178, Aug. 2021.

  - Relevance - **GeN-Foam**. Main reference for one- and two-phase Euler-Euler
    solvers: governing equations, solution algorithms.

.. _RADMAN2021111422-2:

- S. Radman, C. Fiorina, and A. Pautz, "`Development of a novel two-phase flow
  solver for nuclear reactor analysis: Validation against sodium boiling experiments
  <https://doi.org/10.1016/j.nucengdes.2021.111422>`_", Nuclear Engineering and
  Design, vol. 384, p. 111 422, 2021, ISSN: 0029-5493.

  - Relevance - **GeN-Foam**. Important reference for one- and two-phase Euler-Euler
    solvers: treatment of sub-scale structures, interaction models
    (fluid-fluid and fluid-structure), flow regime map modelling, physical models
    (with a focus on sodium).

.. _QIN2023112320:

- H. Qin, C. Fiorina, R. Zhang, S. Radman, D. Zhang, C. Wang, W. Tian, S. Qiu, and G. Su,
  "`Extension of GeN-Foam to modeling of boiling water and validation against the OECD/NRC
  PSBT benchmark <https://doi.org/10.1016/j.nucengdes.2023.112320>`_,"
  Nuclear Engineering and Design, vol. 408, p. 112320, 2023.

  - Relevance - **GeN-Foam**. Reference publication for two-phase flow modeling of water,
    before departure from nucleate boiling.

.. _NIU2024112748:

- Y. Niu, S. Alessandro, C. Fiorina, H. Qin, G. Lazare, Y. Wu, W. Tian, and G. Su,
  "`Extension of GeN-Foam to departure from nucleate boiling prediction and validation against the
  OECD/NRC PSBT benchmark <https://doi.org/10.1016/j.nucengdes.2023.112748>`_,"
  Nuclear Engineering and Design, vol. 416, p. 112748, 2024.

  - Relevance - **GeN-Foam**. Reference publication for two-phase flow modeling of water,
    focusing on prediction of departure from nucleate boiling  and post-CHF flow conditions.

.. _RADMAN2021111422:

- S. Radman, C. Fiorina, K. Mikityuk, and A. Pautz,
  "`Development of a novel two-phase flow solver
  for nuclear reactor analysis: algorithms, verification and implementation in OpenFOAM
  <https://doi.org/10.1016/j.nucengdes.2021.111178`_", Nuclear Engineering and
  Design, vol. 335, 111178, Dec. 2019.

  - Relevance - **GeN-Foam**. Example of use of thermal and pressure baffles.

.. _FIORINA2019376:

- C. Fiorina, "`Impact of the volume heat source on the RANS-based CFD analysis of
  Molten Salt Reactors<https://doi.org/10.1016/j.anucene.2019.06.024>`_",
  Annals of Nuclear Energy, Vol. 134,  pp. 376-382, Dec. 2019.

  - Relevance - **GeN-Foam**. Description of the thermal wall function for
    volumetrically heated fluids.


----------
Neutronics
----------

- C. Fiorina, I. Clifford, M. Aufiero, and K. Mikityuk, "`GeN-Foam: A novel
  OpenFOAM® based multi-physics solver for 2D/3D transient analysis of nuclear
  reactors <https://doi.org/10.1016/j.nucengdes.2015.05.035>`_",
  Nuclear Engineering and Design, vol. 294, pp. 24-37, Dec. 2015.

  - Relevance - **GeN-Foam**. Basic description of diffusion and DNP equations.

.. _FIORINA2016212:

- C. Fiorina, N. Kerkar, K. Mikityuk, P. Rubiolo, and A. Pautz, "`Development and
  verification of the neutron diffusion solver for the GeN-Foam multi-physics
  platform <https://doi.org/10.1016/j.anucene.2016.05.023>`_", Annals of Nuclear
  Energy, vol. 96, pp. 212-222, 2016, ISSN: 0306-4549.

  - Relevance - **GeN-Foam**. Detailed description of diffusion and DNP equations.

  - Warning - Cross-section parametrization. There are no dedicated
    publications but documentation and tutorials  cover this.

.. _FIORINA2017419:

- C. Fiorina, M. Hursin, A. Pautz. "`Extension of the GeN-Foam neutronic solver to SP3
  analysis and application to the CROCUS experimental reactor
  <https://doi.org/10.1016/j.anucene.2016.11.042>`_", Annals of Nuclear Energy, vol. 101,
  pp. 419-428, March 2017.

  - Relevance - **GeN-Foam**. Detailed description of the SP3 model.

  - Warning - Cross-section parametrization. There are no dedicated
    publications but documentation and tutorials cover this.

.. _MATTIOLI:

- A.S. Mattioli, C. Fiorina, S. Lorenzi, A. Cammi. "`Derivation and implementation in
  OpenFOAM of a point-kinetics model for Molten Salt Reactors
  <https://re.public.polimi.it/handle/11311/1225555>`_", ANS Winter Meeting 2021

  - Relevance - **GeN-Foam**. Detailed theory of the PRKE model for MSRs.

.. _RADMAN2022108891:

- S. Radman, C. Fiorina, P. Song, and A. Pautz, "`Development of a point-kinetics
  model in OpenFOAM, integration in GeN-Foam, and validation against FFTF
  experimental data <https://doi.org/10.1016/j.anucene.2021.108891>`_",
  Annals of Nuclear Energy, vol. 168, p. 108 891, 2022, ISSN: 0306-4549.

  - Relevance - **GeN-Foam**. Detailed theory of the PRKE model.

  - Warning - The model has slightly evolved over time. There are no dedicated
    publications but documentation and tutorials should cover its use.

.. _GUILBAUD202310490:

- T. Guilbaud, C. Fiorina, A. Pautz, and F. Carminati, "`Preliminary safety analysis of
  the Transmutex sub-critical reactor using the GeN-Foam multi-physics solver
  <https://doi.org/10.1016/j.pnucene.2023.104980>`_," Progress in Nuclear
  Energy, vol. 168, p. 104980, 2024.

  - Relevance - **GeN-Foam**. Description of the extension of the point-kinetics
    sub-solver for the analysis of sub-critical reactor configurations.


-----------------------------------------
Structural mechanics and fuel performance
-----------------------------------------

.. _BRUNETTO2023112:

- E. L. Brunetto, A. Scolaro, C. Fiorina, and A. Pautz, "`Extension of the OFFBEAT
  fuel performance code to finite strains and validation against LOCA experiments
  <https://doi.org/10.1016/j.nucengdes.2023.112232>`_", Nuclear Engineering and
  Design, vol. 406, pp. 112-232, 2023, ISSN: 0029-5493.

  - Relevance - **OFFBEAT**. Description of models for finite strains.

.. _SCOLARO2022309:

- A. Scolaro, C. Fiorina, I. Clifford, and A. Pautz, "`Development of a semi-implicit contact
  methodology for finite volume stress solvers <https://doi.org/10.1002/nme.6857>`_,"
  International Journal for Numerical Methods in Engineering, vol. 123, no. 2,
  pp. 309-338, 2022.

  - Relevance - **OFFBEAT**. Description of the implicit contact boundary condition.


---------------------
Multiphysics coupling
---------------------

- C. Fiorina, I. Clifford, M. Aufiero, and K. Mikityuk, "`GeN-Foam: A novel
  OpenFOAM® based multi-physics solver for 2D/3D transient analysis of nuclear
  reactors <https://doi.org/10.1016/j.nucengdes.2015.05.035>`_",
  Nuclear Engineering and Design, vol. 294, pp. 24-37, Dec. 2015.

  - Relevance - **GeN-Foam**. Description of algorithms for mesh deformation and
    connection between neutronics and thermal-mechanics.

.. _FIORINA201925:

- C. Fiorina, S. Radman, M.-Z. Koc, and A. Pautz, "`Detailed modelling of the expansion
  reactivity feedback in fast reactors using OpenFOAM
  <https://www.researchgate.net/publication/337739496_DETAILED_MODELLING_OF_THE_EXPANSION_REACTIVITY_FEEDBACK_IN_FAST_REACTORS_USING_OpenFOAM>`_",
  in Proceedings of the International Conference on Mathematics and Computational Methods
  Applied to Nuclear Science and Engineering, M and C, Portland, OR, USA, pp. 25-29, 2019.

  - Relevance - **GeN-Foam**. Detailed example of a case involving the prediction of
    reactivity feedback associated with thermal expansion.


------------
FMI coupling
------------

.. _GUILBAUD2024105022:

- T. Guilbaud, et al. "`Investigating the Functional Mock-up Interface as a Coupling
  Framework for the multi-fidelity analysis of nuclear reactors
  <https://doi.org/10.1016/j.pnucene.2023.105022>`_", Progress in Nuclear Energy
  vol. 169, 105022, April 2024.

  - Relevance - **foamForNuclear**. Description of the implementation and use
    of the FMI interface.


---------------------------
Other relevant publications
---------------------------

- C. Fiorina, I. Clifford, S. Kelm, and S. Lorenzi, "`On the development of multi-
  physics tools for nuclear reactor analysis based on OpenFOAM (R): state of the
  art, lessons learned and perspectives <https://doi.org/10.1016/j.nucengdes.2021.111604>`_",
  Nuclear Engineering And Design, vol. 387, p. 111 604, 2022.

  - Relevance - **GeN-Foam**. Overview on the use of OpenFOAM for nuclear engineering
    applications.

- P. German, et al., "`GeN-ROM—An OpenFOAM®-based multiphysics reduced-order modeling
  framework for the analysis of Molten Salt Reactors
  <https://doi.org/10.1016/j.pnucene.2022.104148>`_", Progress in Nuclear Energy
  Vol. 146, 104148,  April 2022.

  - Relevance - **GeN-Foam**. Library for POD-based reduced-order modeling in GeN-Foam

  - Warning. This is is a separate project, not maintained by foamForNuclear. Its compatibility
    with the current foamForNuclear libraries is not guaranteed.


------------------
Learning resources
------------------

- "`Multi-physics modelling and simulation of nuclear reactors using OpenFOAM
   <https://elearning.iaea.org/m2/course/view.php?id=1286>`_"

  - Relevance - **GeN-Foam, OFFBEAT, OpenFOAM**. Introduction to OpenFOAM,
    GeN-Foam, OFFBEAT.

  - Warning. Requires an IAEA Nucleus account. The course was offered prior
    to foamForNuclear, so prior to a new structure that merges the capabilities
    of GeN-Foam and OFFBEAT, and widely expand the multi-physics capabilities
    of the platform. However, it still provides a well structured
    introduction on how to approach complex OpenFOAM-based tools like foamForNuclear



TO BE ADDED (or not) @Ale

- Investigation on the effect of eccentricity for fuel disc irradiation tests
  A Scolaro, P Van Uffelen, C Fiorina, A Schubert, I Clifford, A Pautz
  Nuclear Engineering and Technology 53 (5), 1602-1611

- Coupling methodology for the multidimensional fuel performance code offbeat and the Monte Carlo neutron transport code SERPENT
  A Scolaro, Y Robert, C Fiorina, I Clifford, A Pautz
  Proceedings of the Global/Top Fuel

- Cladding plasticity modeling with the multidimensional fuel performance code OFFBEAT
  A Scolaro, I Clifford, C Fiorina, A Pautz
  Global/TopFuel

- Simulation of fission gas release in the 3D fuel performance code OFFBEAT
  A Gianesello, A Scolaro, C Fiorina
  IAEA Fast Reactor Conference - FR22

- International fuel performance study of fresh fuel experiments for PCMI effects during RIA experiments
  S Seo, C Folsom, C Jensen, D Kamerman, L Giaccardi, M Cherubini, ...
  Nuclear Engineering and Design 430, 113673

- PuMMA blind benchmark: Performance of high plutonium content MOX fuel under irradiation
  D Jaramillo-Sierra, M Stefanowska-Skrodzka, J Lavarenne, E Deveaux, ...
  Nuclear Engineering and Design 435, 113960

- Application of the OFFBEAT Fuel Performance Code to IAEA Verification Exercises on LMFBR Core Bowing
  A Cornet, E Brunetto, A Scolaro, C Fiorina, A Pautz
  The International Conference on Mathematics and Computational Methods …

- Development of a dynamic-mesh porosity transport model for multi-dimensional fuel performance codes
  EL Brunetto, C Fiorina, A Pautz, S van Til, F Nindiyasari, A Fedorov, ...
  Journal of Nuclear Materials 608, 155717

- Modeling of Zircaloy Oxidation Through Dynamic Mesh Deformation
  A Scolaro, E Brunetto, C Fiorina
  arXiv preprint arXiv:2404.03454

- Implementation of multi-dimensional transport solvers in OFFBEAT
  EL Brunetto, M Reymond, A Scolaro, C Fiorina, A Pautz
  Nufuel 2023: Book of abstracts, 12-12

- Preliminary Extension of OFFBEAT to TRISO Fuel
  F Xiang, A Scolaro, Y Wu, S Qiu, C Fiorina, A Pautz
  Nufuel 2023: Book of abstracts, 6-6

- Fuel pin behaviour under irradiation with high Pu content: Benchmark exercise
  M Stefanowska-Skrodzka, J Lavarenne, E Deveaux, E Brunetto, ...
