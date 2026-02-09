.. _vnv_references:

==============
V&V References
==============


-------------------------------------------
Thermal-hydraulics / heat and mass transfer
-------------------------------------------

- C. Fiorina, I. Clifford, M. Aufiero, and K. Mikityuk, "`GeN-Foam: A novel
  OpenFOAM® based multi-physics solver for 2D/3D transient analysis of nuclear
  reactors <https://doi.org/10.1016/j.nucengdes.2015.05.035>`_",
  Nuclear Engineering and Design, vol. 294, pp. 24-37, Dec. 2015.

  - Relevance - **GeN-Foam**.  Comparison with analytic results for 1-D
    simulations of porous structures. Comparison with analytic results for 1-D simulations
    using a local coordinate system for the porous structure that is rotated compared to the
    global coordinate system. Consistency check based on natural convection cavity flow using
    both the compressible version of the sub-solver and the incompressible version with
    Boussinesq approximation. Verification against analytical solution of the sub-scale
    structure for nuclear fuel.

- S. Radman, C. Fiorina, and A. Pautz, "`Development of a novel two-phase flow solver
  for nuclear reactor analysis: algorithms, verification and implementation in OpenFOAM
  <https://doi.org/10.1016/j.nucengdes.2021.111178>`_", Nuclear Engineering and
  Design, vol. 379, 111178, Aug. 2021.

  - Relevance - **GeN-Foam**. Verification for one- and two-phase Euler-Euler solvers
    using the method of manufactured solutions


- S. Radman, C. Fiorina, and A. Pautz, "`Development of a novel two-phase flow
  solver for nuclear reactor analysis: Validation against sodium boiling experiments
  <https://doi.org/10.1016/j.nucengdes.2021.111422>`_", Nuclear Engineering and
  Design, vol. 384, p. 111 422, 2021, ISSN: 0029-5493.

  - Relevance - **GeN-Foam**. Validation or one- and two-phase Euler-Euler solvers.
    This paper focused on validation effort against steady state
    and transient sodium boiling experiments. The steady-state experiments,
    performed at JRC, ISPRA, consisted in sodium boiling in
    heated tube geometries at constant inlet flow velocities and inlet power
    levels until a quasi-steady state is reached. The results of this validation
    phase were satisfactory. Test L22 at the KNS-37 facility at Kernforschungszentrum
    Karlsruhe was then modelled and reanalyzed to preliminary validate the transient
    capabilities of the solver. Different
    correlations for a variety of quantities where investigated, some of
    which were proven to yield satisfactory results when compared to the
    experimental ones.


- H. Qin, C. Fiorina, R. Zhang, S. Radman, D. Zhang, C. Wang, W. Tian, S. Qiu, and G. Su,
  "`Extension of GeN-Foam to modeling of boiling water and validation against the OECD/NRC
  PSBT benchmark <https://doi.org/10.1016/j.nucengdes.2023.112320>`_,"
  Nuclear Engineering and Design, vol. 408, p. 112320, 2023.

  - Relevance - **GeN-Foam**. Validation against the OECD/NEA PSBT benchmark for water boiling.


- Y. Niu, S. Alessandro, C. Fiorina, H. Qin, G. Lazare, Y. Wu, W. Tian, and G. Su,
  "`Extension of GeN-Foam to departure from nucleate boiling prediction and validation against the
  OECD/NRC PSBT benchmark <https://doi.org/10.1016/j.nucengdes.2023.112748>`_,"
  Nuclear Engineering and Design, vol. 416, p. 112748, 2024.

  - Relevance - **GeN-Foam**. Validation against the OECD/NEA PSBT benchmark for water boiling.

- T. Guilbaud, E. Simonnot, A. Scolaro, and C. Fiorina, "`Full core study of the
  KIWI-B-4E Nuclear Thermal Propulsion system using OpenMC and GeN-Foam
  <https://doi.org/10.1016/j.nucengdes.2024.113639>`_", Nuclear Engineering and
  Design, vol. 429, p. 113 639, 2024, ISSN: 0029-5493.

  - Relevance - **GeN-Foam**. Validation of porous-medium thermal hydraulics and
    the hydrogen thermo-physical properties against data from the  KIWI-B-4E experiment

- C. Fiorina, "`Impact of the volume heat source on the RANS-based CFD analysis of
  Molten Salt Reactors<https://doi.org/10.1016/j.anucene.2019.06.024>`_",
  Annals of Nuclear Energy, Vol. 134,  pp. 376-382, Dec. 2019.

  - Relevance - **GeN-Foam**. Verification and validation of
    the thermal wall function for volumetrically heated fluids.

- Z. Hughes, et al., "`Investigating the accuracy of porous-medium treatments in the analysis of
  nuclear thermal propulsion systems<https://www.ans.org/pubs/proceedings/article-55589/>`_",
  Proceedings of the International Conference on Physics of Reactors (PHYSOR 2024)

  - Relevance - **GeN-Foam**. VErification of porous-medium models against heterogeneous
    models for nuclear thermal propulsion.


----------
Neutronics
----------

.. _FIORINA201515226:

- Carlo Fiorina, Konstantin Mikityuk, "Application of the new GeN-Foam
  multi-physics solver to the European Sodium Fast Reactor and verification
  against available codes", Proceedings of ICAPP 2015, May 03-06, 2015 -
  Nice (France), Paper 15226.

- C. Fiorina, I. Clifford, M. Aufiero, and K. Mikityuk, "`GeN-Foam: A novel
  OpenFOAM® based multi-physics solver for 2D/3D transient analysis of nuclear
  reactors <https://doi.org/10.1016/j.nucengdes.2015.05.035>`_",
  Nuclear Engineering and Design, vol. 294, pp. 24-37, Dec. 2015.

  - Relevance - **GeN-Foam**. Verification of the diffusion solver against Serpent
    for an SFR.

- C. Fiorina, N. Kerkar, K. Mikityuk, P. Rubiolo, and A. Pautz, "`Development and
  verification of the neutron diffusion solver for the GeN-Foam multi-physics
  platform <https://doi.org/10.1016/j.anucene.2016.05.023>`_", Annals of Nuclear
  Energy, vol. 96, pp. 212-222, 2016, ISSN: 0306-4549.

  - Relevance - **GeN-Foam**. Verification of the diffusion solver against
    Monte Carlo results for both a thermal and a fast reactor system.

- C. Fiorina, M. Hursin, A. Pautz. "`Extension of the GeN-Foam neutronic solver to SP3
  analysis and application to the CROCUS experimental reactor
  <https://doi.org/10.1016/j.anucene.2016.11.042>`_", Annals of Nuclear Energy, vol. 101,
  pp. 419-428, March 2017.

  - Relevance - **GeN-Foam**. Verification of the SP3 solver against
    Monte Carlo results for both a PWR and the CROCUS reactor.

- C. Fiorina, S. Radman, M.-Z. Koc, and A. Pautz, "`Detailed modelling of the expansion
  reactivity feedback in fast reactors using OpenFOAM
  <https://www.researchgate.net/publication/337739496_DETAILED_MODELLING_OF_THE_EXPANSION_REACTIVITY_FEEDBACK_IN_FAST_REACTORS_USING_OpenFOAM>`_",
  in Proceedings of the International Conference on Mathematics and Computational Methods
  Applied to Nuclear Science and Engineering, M and C, Portland, OR, USA, pp. 25-29, 2019.

  - Relevance - **GeN-Foam**. Validation of of the Sn solver against Godiva. Verification
    against Serpent-OpenFOAM coupling.

- A.S. Mattioli, C. Fiorina, S. Lorenzi, A. Cammi. "`Derivation and implementation in
  OpenFOAM of a point-kinetics model for Molten Salt Reactors
  <https://re.public.polimi.it/handle/11311/1225555>`_"

  - Relevance - **GeN-Foam**. Verification of the PRKE model for MSRs against diffusion.

- T. Mager, C. Fiorina, M. Hursin, and A. Pautz, "`Validation of the GeN-Foam model of the
  Crocus experimental reactor <https://doi.org/10.13182/xyz-33754>`_,"
  in proceedings of the international conference on mathematics
  and computational methods applied to nuclear science and engineering, M and C 2021, 2021.

  - Relevance - **GeN-Foam**. Verification of the Sn solver against
    Monte Carlo results and experimental data.

- T. Mager, O. Pakari, M. Hursin, T. Ligonnet, A. Pautz, and V. Lamirand, "`Intra-pin
  reaction rate measurements in the Crocus experimental reactor
  <https://doi.org/10.1109/TNS.2024.3370371>`_," IEEE Transactions on Nuclear
  Science, 2024.

  - Relevance - **GeN-Foam**. Verification of the Sn solver against
    Monte Carlo results and experimental data.

- B. Lindley, et al., "`Impact of thermal-hydraulic feedback and differential thermal expansion on European
  sodium-cooled fast reactor core power distribution<https://doi.org/10.1115/1.4056930>`_,"
  Journal of Nuclear Engineering and Radiation Science 9 (3), 031301

  - Relevance - **GeN-Foam**. Verification of GeN-Foam neutronics against other codes.

- E. Nikitin, et al., "`Neutronic modelling of the FFTF control rod worth measurements with diffusion codes
  <https://www.epj-conferences.org/articles/epjconf/pdf/2021/01/epjconf_physor2020_10017.pdf>`_,"
  EPJ Web of Conferences 247, 10017

  - Relevance - **GeN-Foam**. Verification of GeN-Foam neutronics against experimental data and other codes
    for the FFTF.

- T. Mager, et al., "`Comparison of high-fidelity transport code predictions against the NECTAR intra-pin
  reaction rate measurements<https://doi.org/10.1016/j.anucene.2025.111430>`_,"
  Annals of Nuclear Energy 219, 111430

  - Relevance - **GeN-Foam**. Validation of GeN-Foam neutronics against experimental data for the
    CROCUS reactor.


----------------------------------------
Multiphysics coupling and integral tests
----------------------------------------

- C. Fiorina, I. Clifford, M. Aufiero, and K. Mikityuk, "`GeN-Foam: A novel
  OpenFOAM® based multi-physics solver for 2D/3D transient analysis of nuclear
  reactors <https://doi.org/10.1016/j.nucengdes.2015.05.035>`_",
  Nuclear Engineering and Design, vol. 294, pp. 24-37, Dec. 2015.

  - Relevance - **GeN-Foam**. Code-to-code benchmark against TRACE for an SFR core
    (steady state and transient)

- C. Fiorina, K. Mikityuk, "`Application of the new GeN-Foam multi-physics solver to
  the European Sodium Fast Reactor and verification against available codes
  <https://www.researchgate.net/profile/Carlo-Fiorina/publication/287923217_Application_of_the_new_GeN-Foam_multi-physics_solver_to_the_European_Sodium_Fast_Reactor_and_verification_against_available_codes/links/61e1349f8d338833e368eb8e/Application-of-the-new-GeN-Foam-multi-physics-solver-to-the-European-Sodium-Fast-Reactor-and-verification-against-available-codes.pdf>`_",
  Proceedings of ICAPP 2015, May 03-06, 2015.

  - Relevance - **GeN-Foam**. Code-to-code benchmark against TRACE for an SFR core
    (steady state and transient)

- S. Radman, C. Fiorina, P. Song, and A. Pautz, "`Development of a point-kinetics
  model in OpenFOAM, integration in GeN-Foam, and validation against FFTF
  experimental data <https://doi.org/10.1016/j.anucene.2021.108891>`_",
  Annals of Nuclear Energy, vol. 168, p. 108 891, 2022, ISSN: 0306-4549.

  - Relevance - **GeN-Foam**. Validation of GeN-Foam against steady-state and
    transient data from the FFTF

.. _GUILBAUD202313182:

- T. Guilbaud and C. Fiorina, "`Multi-physics modeling of the KIWI-B-4E Nuclear
  Thermal Propulsion system using OpenMC and GeN-Foam: preliminary as-
  sessment against experimental data <https://doi.org/10.13182/T128-42091>`_",
  In: Transactions of the American Nuclear
  Society, Indianapolis, Indiana, USA (June 11-14, 2023)

  - Relevance - **GeN-Foam**. Validation against data from the KIWI-B-4E nuclear
    thermal propulsion reactor.

- T. Guilbaud, E. Simonnot, A. Scolaro, and C. Fiorina, "`Full core study of the
  KIWI-B-4E Nuclear Thermal Propulsion system using OpenMC and GeN-Foam
  <https://doi.org/10.1016/j.nucengdes.2024.113639>`_", Nuclear Engineering and
  Design, vol. 429, p. 113 639, 2024, ISSN: 0029-5493.

  - Relevance - **GeN-Foam**. Validation against data from the KIWI-B-4E nuclear
    thermal propulsion reactor.

- R. G. De Oliveira and K. Mikityuk, "`Analytical solutions to a coupled fluid dynamics and
  neutron transport problem with application to GeN-Foam verification
  <https://doi.org/10.1016/j.anucene.2018.07.036>`_," Annals of Nuclear Energy,
  vol. 121, pp. 446-451, 2018.

  - Relevance - **GeN-Foam**. Verification against the CNRS benchmark for molten
    salt reactors.

- M. Tiberga, R. G. G. de Oliveira, E. Cervi, J. A. Blanco, S. Lorenzi,
  M. Aufiero, D. Lathouwers, and P. Rubiolo, "`Results from a multi-physics
  numerical benchmark for codes dedicated to molten salt fast reactors
  <https://doi.org/10.1016/j.anucene.2020.107428>`_," Annals of Nuclear Energy,
  vol. 142, p. 107428, 2020.

  - Relevance - **GeN-Foam**. Verification against the CNRS benchmark for molten
    salt reactors.


- T. Guilbaud, C. Fiorina, A. Pautz, and F. Carminati, "`Preliminary safety analysis of
  the Transmutex sub-critical reactor using the GeN-Foam multi-physics solver
  <https://doi.org/10.1016/j.pnucene.2023.104980>`_," Progress in Nuclear
  Energy, vol. 168, p. 104980, 2024.

  - Relevance - **GeN-Foam**. Verification of the point reactor kinetics model against
    analytic solutions

.. _HABTEMARIAM2024110237:

- N. Habtemariam, C. Fiorina, S. Lorenzi, A. Cammi. "`On the need for multi-dimensional models
  for the safety analysis of (fast-spectrum) Molten Salt Reactors
  <https://doi.org/10.1016/j.anucene.2023.110237>`_", Annals of Nuclear Energy 197, 110237,
  March 2024.

  - Relevance - **GeN-Foam**. Comparison among GeN-Foam models based on different modeling choices.


------
Others
------

- C. Fiorina, "`The GeN-Foam Multiphysics Solver: A Status Update
  <http://dx.doi.org/10.13182/PHYSOR22-37774>`_", in International Conference on
  Physics of Reactors 2022 (PHYSOR 2022), Pittsburgh, PA, Aug. 2022.

  - Relevance - **GeN-Foam**. Summary of the V&V status in 2022.


--------------------------------------------
Mechanics and fuel behavior
--------------------------------------------

- A. Scolaro, I. Clifford, C. Fiorina, and A. Pautz, "`First steps towards the
  development of a 3D nuclear fuel behavior solver with OpenFOAM <https://doi.org/10.1115/ICONE26-82381>`_",
  in Proceedings of the International Conference on Nuclear Engineering (ICONE),
  vol. 51456, paper V003T02A051, ASME, 2018.

  - Relevance - **OFFBEAT**. First thermo-mechanical verification steps for a novel
    multi-dimensional OpenFOAM-based fuel performance solver.


- A. Scolaro, I. Clifford, C. Fiorina, and A. Pautz, "`Multi-dimensional creep analysis
  using the novel OFFBEAT fuel performance code 
  <https://www.iaea.org/publications/14814/progress-on-pellet-cladding-interaction-and-stress-corrosion-cracking>`_",
  in IAEA Technical Meeting on Progress on Pellet Cladding Interaction and Stress Corrosion Cracking,
  2019.

  - Relevance - **OFFBEAT**. First implementation and verification of creep modelling in a multi-dimensional
    finite-volume fuel performance context.


- A. Scolaro, Y. Robert, C. Fiorina, I. Clifford, and A. Pautz, "Coupling methodology for the
  multidimensional fuel performance code OFFBEAT and the Monte Carlo neutron transport code SERPENT",
  in Proceedings of the Global/Top Fuel, 2019.

  - Relevance - **OFFBEAT**. Coupling methodology between Monte Carlo neutronics and fuel
    behaviour, using the Serpent Multi-physics interface. Verification against results obtained with TURBNP module.


- A. Scolaro, I. Clifford, C. Fiorina, and A. Pautz, "Cladding plasticity modeling with the
  multidimensional fuel performance code OFFBEAT",
  in Proceedings of Global/TopFuel, 2019.

  - Relevance - **OFFBEAT**. Verification of elasto-plastic constitutive modelling in the
    finite-volume solid mechanics framework and assessment of stress/strain localization in PCMI.


- A. Scolaro, I. Clifford, C. Fiorina, and A. Pautz, "`The OFFBEAT multi-dimensional fuel behavior solver
  <https://doi.org/10.1016/j.nucengdes.2019.110416>`_",
  Nuclear Engineering and Design, vol. 358, 110416, 2020.

  - Relevance - **OFFBEAT**. Core reference for the solver architecture and main modelling
    capabilities; verification on representative fuel-rod thermo-mechanical problems and
    initial validation against IFA432 and IFA562 rods.


- A. Scolaro, P. Van Uffelen, C. Fiorina, A. Schubert, I. Clifford, and A. Pautz, "`Investigation on the effect
  of eccentricity for fuel disc irradiation tests
  <https://doi.org/10.1016/j.net.2020.11.003>`_",
  Nuclear Engineering and Technology, vol. 53, no. 5, pp. 1602--1611, 2021.

  - Relevance - **OFFBEAT**. Validation study against data from the HRBP disc irradiation tests. 
    The parametric analysis of the differences agaisnt the high-temperature cases support the need for multi-dimensional
    modelling in disc irradiation tests.


- A. Scolaro, C. Fiorina, I. Clifford, and A. Pautz, "`Development of a semi-implicit contact methodology for
  finite volume stress solvers <https://doi.org/10.1002/nme.6857>`_",
  International Journal for Numerical Methods in Engineering, vol. 123, no. 2,
  pp. 309--338, 2022.

  - Relevance - **OFFBEAT**. Verification and code-to-code comparison of the semi-implicit OFFBEAT
    contact algorithm for solid mechanics, including patch tests,
    Hertzian contact, and punch tests.


- E. L. Brunetto, A. Scolaro, C. Fiorina, and A. Pautz, "`Extension of the OFFBEAT fuel performance code to finite strains and validation against LOCA experiments
  <https://doi.org/10.1016/j.nucengdes.2023.112232>`_",
  Nuclear Engineering and Design, vol. 406, 112232, 2023.

  - Relevance - **OFFBEAT**. Verification of finite-strain large-deformation
    thermo-mechanics and validation against integral LOCA experiments (IFA650.2).


- G. Zullo, D. Pizzocri, A. Scolaro, P. Van Uffelen, F. Feria, L. E. Herranz, and L. Luzzi,
  "`Integral-scale validation of the SCIANTIX code for Light Water Reactor fuel rods
  <https://doi.org/10.1016/j.jnucmat.2024.155305>`_",
  *Journal of Nuclear Materials*, vol. 601, 155305, 2024.

  - Relevance - **OFFBEAT**. Comprehensive validation activity for thermal related quantities (fuel centerline
    and fission gas release fraction) against several rods from the IFPE database (IFA432, REGATE, CONTACT, SuperRamp, Risoe3, etc.).


- S. Seo, C. Folsom, C. Jensen, D. Kamerman, L. Giaccardi, M. Cherubini, P. Suk,
  M. Sevecek, J. Sercombe, I. Guenot-Delahaie, et al.,
  "`International fuel performance study of fresh fuel experiments for PCMI effects during RIA experiments
  <https://doi.org/10.1016/j.nucengdes.2024.113673>`_",
  *Nuclear Engineering and Design*, vol. 430, 113673, 2024.

  - Relevance - **OFFBEAT**. Blind international benchmark for RIA PCMI based on
    TREAT and NSRR experiments; cross-code comparison of cladding deformation
    and failure sensitivity to pulse width and deposited energy.


- L. Verma, I. Clifford, P. Konarski, A. Scolaro, and H. Ferroukhi, "`OFFBEAT V\&V studies for REBEKA tests on
  cladding ballooning and burst during LOCA conditions <https://doi.org/10.1016/j.anucene.2024.110773>`_",
  Annals of Nuclear Energy, vol. 208, 110773, 2024.

  - Relevance - **OFFBEAT**. Validation against REBEKA integral LOCA tests for cladding ballooning/burst.


- Q. Faure, G. Delipei, A. Scolaro, M. Avramova, and K. Ivanov, "`Fuel performance code to code comparative analysis
  for the OECD/NEA MPCMIV benchmark
  <https://doi.org/10.1016/j.nucengdes.2024.113685>`_",
  Nuclear Engineering and Design, vol. 430, 113685, 2024.

  - Relevance - **OFFBEAT**. Code-to-code and experimental benchmark for PCMI during
    base irradiation and fast ramp transient; assessment of cladding deformation
    and sensitivity to multi-physics-derived power histories.


- M. Reymond, J. Sercombe, and A. Scolaro, "`Investigation of the PCMI failure of pre-hydrided Zy-4 cladding during
  Reactivity Initiated Accidents with ALCYONE and OFFBEAT fuel performance codes
  <https://doi.org/10.1016/j.nucengdes.2024.113430>`_",
  *Nuclear Engineering and Design*, vol. 427, 113430, 2024.

  - Relevance - **OFFBEAT**. Code-to-code and experimental validation of RIA PCMI
    failure for hydrided cladding; assessment of anisotropic viscoplastic response
    and strain-based failure criteria against ALCYONE and NSRR data.


- L. Verma, I. Clifford, P. Konarski, A. Scolaro, and H. Ferroukhi, "`Analysing hydrogen behaviour in liner claddings
  using OFFBEAT fuel performance code <https://doi.org/10.1016/j.anucene.2025.111559>`_",
  Annals of Nuclear Energy, vol. 222, 111559, 2025.

  - Relevance - **OFFBEAT**. Validation of hydrogen diffusion and hydride-driven
    redistribution in liner claddings against out-of-pile experiments,
    including migration toward the substrate–liner interface.


- E. L. Brunetto, C. Fiorina, A. Pautz, S. van Til, F. Nindiyasari, A. Fedorov, and A. Scolaro,
  "`Development of a dynamic-mesh porosity transport model for multi-dimensional fuel performance codes
  <https://doi.org/10.1016/j.jnucmat.2025.155717>`_",
  Journal of Nuclear Materials, vol. 608, 155717, 2025.

  - Relevance - **OFFBEAT**. Verification and experimental validation of dynamic-mesh
    porosity migration and central hole evolution, including off-centred cavity
    formation and comparison with post-irradiation data.


- D. Jaramillo-Sierra, M. Stefanowska-Skrodzka, J. Lavarenne, E. Deveaux, E. Brunetto,
  V. Matocha, A. Magni, K. Sturm, K. Mikityuk, Y. Wang, et al.,
  "`PuMMA blind benchmark: Performance of high plutonium content MOX fuel under irradiation
  <https://doi.org/10.1016/j.nucengdes.2025.113960>`_",
  *Nuclear Engineering and Design*, vol. 435, 113960, 2025.

  - Relevance - **OFFBEAT**. International blind benchmark for high-Pu MOX fuel;
    cross-code assessment of irradiation behaviour and model uncertainty prior
    to comparison with PIE data.