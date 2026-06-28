.. _usersguide_tutorials:

---------
Tutorials
---------

A set of tutorials is distributed to provide a relatively comprehensive overview
of the capabilities of FFN. In each tutorial:

- a ``README.md`` file provides a general description of the tutorials
- the input files (dictionaries) are normally commented to help identify the role of non-obvious parameters
- an ``Allrun`` bash script is provided to run the tutorial. The steps in the ``Allrun`` script can also help you understand how to run other simulations
- for computationally intensive cases, an ``Allrun_parallel`` bash script is provided to run the tutorial using multiple cores
- for complex cases, an ``Allclean`` script is provided to clean up a case after running it and before starting another simulation


Cases
=====

The Tutorials folder is divided into four types of tutorials:

- `guidedCases <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/guidedCases?ref_type=heads>`_: simple cases to teach the user how to use FFN
- `reactorCases <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases?ref_type=heads>`_: complex cases representing partial or complete nuclear reactors
- `featureCases <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases?ref_type=heads>`_: tests and validation cases for specific features
- `fmuCases <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/fmuCases?ref_type=heads>`_: examples of use of the FMI standard


Regression test suite
=====================

The `Alltest <https://gitlab.com/foamForNuclear/foamForNuclear/-/blob/master/tutorials/Alltest?ref_type=heads>`_ script can be used during development to ensure that
FFN behaves as expected. This script runs all the tutorials contained in this
folder. The process is extensive and may require a day of computation,
depending on the available hardware.

The results of the test are summarized in the
`regressionResults.md <https://gitlab.com/foamForNuclear/foamForNuclear/-/blob/master/tutorials/regressionResults.md?ref_type=heads>`_ file. In ``--fast`` mode, the
results are summarized in
`regressionResults-fast.md <https://gitlab.com/foamForNuclear/foamForNuclear/-/blob/master/tutorials/regressionResults-fast.md?ref_type=heads>`_.

.. code :: bash

    ./Alltest

    # Parallelize the cases
    ./Alltest -j<num_jobs>

    # To add fmi cases to be run, works in parallel
    ./Alltest --fmi

    # To run only fast cases, works in parallel (< 1 min per case)
    ./Alltest --fast

    # To run with the minimum information printed on the terminal
    ./Alltest --quiet



Tutorials
=========

.. raw:: html

   <input type="text" class="tableFilter" placeholder="Search tutorials...">

.. list-table:: Tutorials
    :widths: 50 50 20
    :header-rows: 1
    :class: wy-table-responsive filterable-table

    * - Tutorial title
      - Tags and categories
      - Link to description
    * - 1D Critical Heat Flux
      - | |tutorials/featureCases/1D_CHF-ThermalHydraulics-twoPhase-badge|
      - :doc:`Link <tutorials/featureCases/1D_CHF>`
    * - 1D Heat Exchanger
      - | |tutorials/featureCases/1D_HX-ThermalHydraulics-onePhase-badge|
        | |tutorials/featureCases/1D_HX-ThermalHydraulics-twoPhase-badge|
      - :doc:`Link <tutorials/featureCases/1D_HX>`
    * - 1D PSBT SC
      - | |tutorials/featureCases/1D_PSBT_SC-ThermalHydraulics-twoPhase-badge|
      - :doc:`Link <tutorials/featureCases/1D_PSBT_SC>`
    * - 1D Boiling
      - | |tutorials/featureCases/1D_boiling-ThermalHydraulics-twoPhase-badge|
      - :doc:`Link <tutorials/featureCases/1D_boiling>`
    * - Hron–Turek FSI3 Benchmark
      - | |tutorials/featureCases/2D_HronTurekFSI3-ThermalHydraulics-rhoPimpleFoam-badge|
        | |tutorials/featureCases/2D_HronTurekFSI3-ThermalMechanics-extendedThermoMechanics-badge|
        | |tutorials/featureCases/2D_HronTurekFSI3-Multiphysics-fsiLoop-badge|
      - :doc:`Link <tutorials/featureCases/2D_HronTurekFSI3>`
    * - 2D KNS37-L22
      - | |tutorials/featureCases/2D_KNS37-L22-ThermalHydraulics-twoPhase-badge|
      - :doc:`Link <tutorials/featureCases/2D_KNS37-L22>`
    * - 2D Cavity Boussinesq
      - | |tutorials/featureCases/2D_cavityBoussinesq-ThermalHydraulics-onePhase-badge|
      - :doc:`Link <tutorials/featureCases/2D_cavityBoussinesq>`
    * - 2D external source diffusion
      - | |tutorials/featureCases/2D_externalSourceDiffusion-Neutronics-diffusion-badge|
      - :doc:`Link <tutorials/featureCases/2D_externalSourceDiffusion>`
    * - Vynnycky Conjugate Heat Transfer (CHT) Benchmark
      - | |tutorials/featureCases/2D_flowOverHeatedPlate-ThermalHydraulics-onePhase-badge|
        | |tutorials/featureCases/2D_flowOverHeatedPlate-ThermalMechanics-extendedThermoMechanics-badge|
        | |tutorials/featureCases/2D_flowOverHeatedPlate-Multiphysics-chtLoop-badge|
      - :doc:`Link <tutorials/featureCases/2D_flowOverHeatedPlate>`
    * - Testing for physics coupling
      - | |tutorials/featureCases/2D_fullCoupling-ThermalHydraulics-onePhase-badge|
        | |tutorials/featureCases/2D_fullCoupling-Neutronics-diffusion-badge|
        | |tutorials/featureCases/2D_fullCoupling-ThermalMechanics-legacyThermoMechanics-badge|
        | |tutorials/featureCases/2D_fullCoupling-Multiphysics-looseCoupling-badge|
      - :doc:`Link <tutorials/featureCases/2D_fullCoupling>`
    * - 2D 1-phase & Point-Kinetics Coupling
      - | |tutorials/featureCases/2D_onePhaseAndPointKineticsCoupling-ThermalHydraulics-onePhase-badge|
        | |tutorials/featureCases/2D_onePhaseAndPointKineticsCoupling-Neutronics-pointKinetics-badge|
        | |tutorials/featureCases/2D_onePhaseAndPointKineticsCoupling-Multiphysics-looseCoupling-badge|
      - :doc:`Link <tutorials/featureCases/2D_onePhaseAndPointKineticsCoupling>`
    * - Subcritical Point-Kinetics Solver
      - | |tutorials/featureCases/2D_onePhaseAndSubcriticalPointKineticsCoupling-ThermalHydraulics-onePhase-badge|
        | |tutorials/featureCases/2D_onePhaseAndSubcriticalPointKineticsCoupling-Neutronics-pointKinetics-badge|
        | |tutorials/featureCases/2D_onePhaseAndSubcriticalPointKineticsCoupling-Multiphysics-looseCoupling-badge|
      - :doc:`Link <tutorials/featureCases/2D_onePhaseAndSubcriticalPointKineticsCoupling>`
    * - 2D Void Motion No Phase Change
      - | |tutorials/featureCases/2D_voidMotionNoPhaseChange-ThermalHydraulics-twoPhase-badge|
      - :doc:`Link <tutorials/featureCases/2D_voidMotionNoPhaseChange>`
    * - 2D-wedge ALFRED model with FMU balance of plant
      - | |tutorials/fmuCases/2D_LFRpowerPlant-ThermalHydraulics-onePhase-badge|
        | |tutorials/fmuCases/2D_LFRpowerPlant-Neutronics-pointKinetics-badge|
        | |tutorials/fmuCases/2D_LFRpowerPlant-Multiphysics-looseCoupling-badge|
        | |tutorials/fmuCases/2D_LFRpowerPlant-Multicode-FMI-badge|
      - :doc:`Link <tutorials/fmuCases/2D_LFRpowerPlant>`
    * - 2D Point-Kinetics coupling with FMU
      - | |tutorials/fmuCases/2D_PKCoupleFMI-ThermalHydraulics-onePhase-badge|
        | |tutorials/fmuCases/2D_PKCoupleFMI-Neutronics-pointKinetics-badge|
        | |tutorials/fmuCases/2D_PKCoupleFMI-Multiphysics-looseCoupling-badge|
        | |tutorials/fmuCases/2D_PKCoupleFMI-Multicode-FMI-badge|
      - :doc:`Link <tutorials/fmuCases/2D_PKCoupleFMI>`
    * - Power-temperature-momentum controller
      - | |tutorials/fmuCases/powerTemperatureMomentumControl-ThermalHydraulics-onePhase-badge|
        | |tutorials/fmuCases/powerTemperatureMomentumControl-Multicode-FMI-badge|
      - :doc:`Link <tutorials/fmuCases/powerTemperatureMomentumControl>`
    * - Neutronics Diffusion: Slab Reactor
      - | |tutorials/guidedCases/1_reactorSlab_1D_1Gr_neutronicDiffusion-Neutronics-diffusion-badge|
      - :doc:`Link <tutorials/guidedCases/1_reactorSlab_1D_1Gr_neutronicDiffusion>`
    * - Neutronics Diffusion: Slab Reactor with a reflector
      - | |tutorials/guidedCases/2_reactorSlabReflected_1D_1Gr_neutronicDiffusion-Neutronics-diffusion-badge|
      - :doc:`Link <tutorials/guidedCases/2_reactorSlabReflected_1D_1Gr_neutronicDiffusion>`
    * - ESFR Workshop
      - | 
      - :doc:`Link <tutorials/guidedCases/ESFRWorkshop>`
    * - Test Cases and Tutorials
      - | 
      - :doc:`Link <tutorials/offbeat/testCases>`
    * - 1D MSR Point-Kinetics
      - | |tutorials/reactorCases/1D_MSR_pointKinetics-ThermalHydraulics-onePhase-badge|
        | |tutorials/reactorCases/1D_MSR_pointKinetics-Neutronics-pointKinetics-badge|
        | |tutorials/reactorCases/1D_MSR_pointKinetics-Multiphysics-looseCoupling-badge|
      - :doc:`Link <tutorials/reactorCases/1D_MSR_pointKinetics>`
    * - 1D thermal MSR Point-Kinetics
      - | |tutorials/reactorCases/1D_thermalMSR_pointKinetics-ThermalHydraulics-onePhase-badge|
        | |tutorials/reactorCases/1D_thermalMSR_pointKinetics-Neutronics-pointKinetics-badge|
        | |tutorials/reactorCases/1D_thermalMSR_pointKinetics-Multiphysics-looseCoupling-badge|
      - :doc:`Link <tutorials/reactorCases/1D_thermalMSR_pointKinetics>`
    * - Fast Flux Test Facility — LOFWOS‑13 Test
      - | 
      - :doc:`Link <tutorials/reactorCases/2D_FFTF>`
    * - 2D Molten Salt Fast Reactor
      - | |tutorials/reactorCases/2D_MSFR-ThermalHydraulics-onePhase-badge|
        | |tutorials/reactorCases/2D_MSFR-Neutronics-diffusion-badge|
        | |tutorials/reactorCases/2D_MSFR-Multiphysics-picardLoop-badge|
      - :doc:`Link <tutorials/reactorCases/2D_MSFR>`
    * - 2D MSFR ULOF - Diffusion and Point-Kinetics
      - | |tutorials/reactorCases/2D_MSFR_ULOF_diffusion_and_pk-ThermalHydraulics-onePhase-badge|
        | |tutorials/reactorCases/2D_MSFR_ULOF_diffusion_and_pk-Neutronics-diffusion-badge|
        | |tutorials/reactorCases/2D_MSFR_ULOF_diffusion_and_pk-Multiphysics-picardLoop-badge|
      - :doc:`Link <tutorials/reactorCases/2D_MSFR_ULOF_diffusion_and_pk>`
    * - Pebble Bed Gas-Cooled Reactor Model (HTR-10 Design)
      - | |tutorials/reactorCases/3D_HTR-10-ThermalHydraulics-onePhase-badge|
      - :doc:`Link <tutorials/reactorCases/3D_HTR-10>`
    * - KIWI-B-4E NTP Full Core
      - | |tutorials/reactorCases/3D_KIWI-B-4E-ThermalHydraulics-onePhase-badge|
        | |tutorials/reactorCases/3D_KIWI-B-4E-Neutronics-SP3-badge|
        | |tutorials/reactorCases/3D_KIWI-B-4E-Multiphysics-looseCoupling-badge|
      - :doc:`Link <tutorials/reactorCases/3D_KIWI-B-4E>`
    * - KIWI-B-4E NTP Fuel Assembly
      - | |tutorials/reactorCases/3D_NTPfuelAssembly-ThermalHydraulics-onePhase-badge|
        | |tutorials/reactorCases/3D_NTPfuelAssembly-Neutronics-diffusion-badge|
        | |tutorials/reactorCases/3D_NTPfuelAssembly-Multiphysics-looseCoupling-badge|
      - :doc:`Link <tutorials/reactorCases/3D_NTPfuelAssembly>`
    * - 3D Small ESFR
      - | |tutorials/reactorCases/3D_SmallESFR-ThermalHydraulics-onePhase-badge|
        | |tutorials/reactorCases/3D_SmallESFR-Neutronics-diffusion-badge|
        | |tutorials/reactorCases/3D_SmallESFR-ThermalMechanics-extendedThermoMechanics-badge|
        | |tutorials/reactorCases/3D_SmallESFR-Multiphysics-looseCoupling-badge|
      - :doc:`Link <tutorials/reactorCases/3D_SmallESFR>`
    * - 3D gFHR Pebble Bed Reactor
      - | |tutorials/reactorCases/3D_gFHR-ThermalHydraulics-onePhase-badge|
      - :doc:`Link <tutorials/reactorCases/3D_gFHR>`
    * - Particle-bed Nuclear Thermal Propulsion Reactor
      - | |tutorials/reactorCases/3D_genericParticleBedNTP-ThermalHydraulics-onePhase-badge|
        | |tutorials/reactorCases/3D_genericParticleBedNTP-ThermalMechanics-extendedThermoMechanics-badge|
        | |tutorials/reactorCases/3D_genericParticleBedNTP-Multiphysics-looseCoupling-badge|
      - :doc:`Link <tutorials/reactorCases/3D_genericParticleBedNTP>`
    * - Godiva sphere using discrete ordinate SN
      - | |tutorials/reactorCases/Godiva_SN-Neutronics-SN-badge|
      - :doc:`Link <tutorials/reactorCases/Godiva_SN>`


.. |tutorials/featureCases/1D_CHF-ThermalHydraulics-twoPhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-twoPhase-blue

.. |tutorials/featureCases/1D_HX-ThermalHydraulics-onePhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-onePhase-blue

.. |tutorials/featureCases/1D_HX-ThermalHydraulics-twoPhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-twoPhase-blue

.. |tutorials/featureCases/1D_PSBT_SC-ThermalHydraulics-twoPhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-twoPhase-blue

.. |tutorials/featureCases/1D_boiling-ThermalHydraulics-twoPhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-twoPhase-blue

.. |tutorials/featureCases/2D_HronTurekFSI3-ThermalHydraulics-rhoPimpleFoam-badge| image:: https://img.shields.io/badge/ThermalHydraulics-rhoPimpleFoam-blue

.. |tutorials/featureCases/2D_HronTurekFSI3-ThermalMechanics-extendedThermoMechanics-badge| image:: https://img.shields.io/badge/ThermalMechanics-extendedThermoMechanics-blue

.. |tutorials/featureCases/2D_HronTurekFSI3-Multiphysics-fsiLoop-badge| image:: https://img.shields.io/badge/Multiphysics-fsiLoop-orange

.. |tutorials/featureCases/2D_KNS37-L22-ThermalHydraulics-twoPhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-twoPhase-blue

.. |tutorials/featureCases/2D_cavityBoussinesq-ThermalHydraulics-onePhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-onePhase-blue

.. |tutorials/featureCases/2D_externalSourceDiffusion-Neutronics-diffusion-badge| image:: https://img.shields.io/badge/Neutronics-diffusion-blue

.. |tutorials/featureCases/2D_flowOverHeatedPlate-ThermalHydraulics-onePhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-onePhase-blue

.. |tutorials/featureCases/2D_flowOverHeatedPlate-ThermalMechanics-extendedThermoMechanics-badge| image:: https://img.shields.io/badge/ThermalMechanics-extendedThermoMechanics-blue

.. |tutorials/featureCases/2D_flowOverHeatedPlate-Multiphysics-chtLoop-badge| image:: https://img.shields.io/badge/Multiphysics-chtLoop-orange

.. |tutorials/featureCases/2D_fullCoupling-ThermalHydraulics-onePhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-onePhase-blue

.. |tutorials/featureCases/2D_fullCoupling-Neutronics-diffusion-badge| image:: https://img.shields.io/badge/Neutronics-diffusion-blue

.. |tutorials/featureCases/2D_fullCoupling-ThermalMechanics-legacyThermoMechanics-badge| image:: https://img.shields.io/badge/ThermalMechanics-legacyThermoMechanics-blue

.. |tutorials/featureCases/2D_fullCoupling-Multiphysics-looseCoupling-badge| image:: https://img.shields.io/badge/Multiphysics-looseCoupling-orange

.. |tutorials/featureCases/2D_onePhaseAndPointKineticsCoupling-ThermalHydraulics-onePhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-onePhase-blue

.. |tutorials/featureCases/2D_onePhaseAndPointKineticsCoupling-Neutronics-pointKinetics-badge| image:: https://img.shields.io/badge/Neutronics-pointKinetics-blue

.. |tutorials/featureCases/2D_onePhaseAndPointKineticsCoupling-Multiphysics-looseCoupling-badge| image:: https://img.shields.io/badge/Multiphysics-looseCoupling-orange

.. |tutorials/featureCases/2D_onePhaseAndSubcriticalPointKineticsCoupling-ThermalHydraulics-onePhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-onePhase-blue

.. |tutorials/featureCases/2D_onePhaseAndSubcriticalPointKineticsCoupling-Neutronics-pointKinetics-badge| image:: https://img.shields.io/badge/Neutronics-pointKinetics-blue

.. |tutorials/featureCases/2D_onePhaseAndSubcriticalPointKineticsCoupling-Multiphysics-looseCoupling-badge| image:: https://img.shields.io/badge/Multiphysics-looseCoupling-orange

.. |tutorials/featureCases/2D_voidMotionNoPhaseChange-ThermalHydraulics-twoPhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-twoPhase-blue

.. |tutorials/fmuCases/2D_LFRpowerPlant-ThermalHydraulics-onePhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-onePhase-blue

.. |tutorials/fmuCases/2D_LFRpowerPlant-Neutronics-pointKinetics-badge| image:: https://img.shields.io/badge/Neutronics-pointKinetics-blue

.. |tutorials/fmuCases/2D_LFRpowerPlant-Multiphysics-looseCoupling-badge| image:: https://img.shields.io/badge/Multiphysics-looseCoupling-orange

.. |tutorials/fmuCases/2D_LFRpowerPlant-Multicode-FMI-badge| image:: https://img.shields.io/badge/Multicode-FMI-red

.. |tutorials/fmuCases/2D_PKCoupleFMI-ThermalHydraulics-onePhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-onePhase-blue

.. |tutorials/fmuCases/2D_PKCoupleFMI-Neutronics-pointKinetics-badge| image:: https://img.shields.io/badge/Neutronics-pointKinetics-blue

.. |tutorials/fmuCases/2D_PKCoupleFMI-Multiphysics-looseCoupling-badge| image:: https://img.shields.io/badge/Multiphysics-looseCoupling-orange

.. |tutorials/fmuCases/2D_PKCoupleFMI-Multicode-FMI-badge| image:: https://img.shields.io/badge/Multicode-FMI-red

.. |tutorials/fmuCases/powerTemperatureMomentumControl-ThermalHydraulics-onePhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-onePhase-blue

.. |tutorials/fmuCases/powerTemperatureMomentumControl-Multicode-FMI-badge| image:: https://img.shields.io/badge/Multicode-FMI-red

.. |tutorials/guidedCases/1_reactorSlab_1D_1Gr_neutronicDiffusion-Neutronics-diffusion-badge| image:: https://img.shields.io/badge/Neutronics-diffusion-blue

.. |tutorials/guidedCases/2_reactorSlabReflected_1D_1Gr_neutronicDiffusion-Neutronics-diffusion-badge| image:: https://img.shields.io/badge/Neutronics-diffusion-blue





.. |tutorials/reactorCases/1D_MSR_pointKinetics-ThermalHydraulics-onePhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-onePhase-blue

.. |tutorials/reactorCases/1D_MSR_pointKinetics-Neutronics-pointKinetics-badge| image:: https://img.shields.io/badge/Neutronics-pointKinetics-blue

.. |tutorials/reactorCases/1D_MSR_pointKinetics-Multiphysics-looseCoupling-badge| image:: https://img.shields.io/badge/Multiphysics-looseCoupling-orange

.. |tutorials/reactorCases/1D_thermalMSR_pointKinetics-ThermalHydraulics-onePhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-onePhase-blue

.. |tutorials/reactorCases/1D_thermalMSR_pointKinetics-Neutronics-pointKinetics-badge| image:: https://img.shields.io/badge/Neutronics-pointKinetics-blue

.. |tutorials/reactorCases/1D_thermalMSR_pointKinetics-Multiphysics-looseCoupling-badge| image:: https://img.shields.io/badge/Multiphysics-looseCoupling-orange



.. |tutorials/reactorCases/2D_MSFR-ThermalHydraulics-onePhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-onePhase-blue

.. |tutorials/reactorCases/2D_MSFR-Neutronics-diffusion-badge| image:: https://img.shields.io/badge/Neutronics-diffusion-blue

.. |tutorials/reactorCases/2D_MSFR-Multiphysics-picardLoop-badge| image:: https://img.shields.io/badge/Multiphysics-picardLoop-orange

.. |tutorials/reactorCases/2D_MSFR_ULOF_diffusion_and_pk-ThermalHydraulics-onePhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-onePhase-blue

.. |tutorials/reactorCases/2D_MSFR_ULOF_diffusion_and_pk-Neutronics-diffusion-badge| image:: https://img.shields.io/badge/Neutronics-diffusion-blue

.. |tutorials/reactorCases/2D_MSFR_ULOF_diffusion_and_pk-Multiphysics-picardLoop-badge| image:: https://img.shields.io/badge/Multiphysics-picardLoop-orange

.. |tutorials/reactorCases/3D_HTR-10-ThermalHydraulics-onePhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-onePhase-blue

.. |tutorials/reactorCases/3D_KIWI-B-4E-ThermalHydraulics-onePhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-onePhase-blue

.. |tutorials/reactorCases/3D_KIWI-B-4E-Neutronics-SP3-badge| image:: https://img.shields.io/badge/Neutronics-SP3-blue

.. |tutorials/reactorCases/3D_KIWI-B-4E-Multiphysics-looseCoupling-badge| image:: https://img.shields.io/badge/Multiphysics-looseCoupling-orange

.. |tutorials/reactorCases/3D_NTPfuelAssembly-ThermalHydraulics-onePhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-onePhase-blue

.. |tutorials/reactorCases/3D_NTPfuelAssembly-Neutronics-diffusion-badge| image:: https://img.shields.io/badge/Neutronics-diffusion-blue

.. |tutorials/reactorCases/3D_NTPfuelAssembly-Multiphysics-looseCoupling-badge| image:: https://img.shields.io/badge/Multiphysics-looseCoupling-orange

.. |tutorials/reactorCases/3D_SmallESFR-ThermalHydraulics-onePhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-onePhase-blue

.. |tutorials/reactorCases/3D_SmallESFR-Neutronics-diffusion-badge| image:: https://img.shields.io/badge/Neutronics-diffusion-blue

.. |tutorials/reactorCases/3D_SmallESFR-ThermalMechanics-extendedThermoMechanics-badge| image:: https://img.shields.io/badge/ThermalMechanics-extendedThermoMechanics-blue

.. |tutorials/reactorCases/3D_SmallESFR-Multiphysics-looseCoupling-badge| image:: https://img.shields.io/badge/Multiphysics-looseCoupling-orange

.. |tutorials/reactorCases/3D_gFHR-ThermalHydraulics-onePhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-onePhase-blue

.. |tutorials/reactorCases/3D_genericParticleBedNTP-ThermalHydraulics-onePhase-badge| image:: https://img.shields.io/badge/ThermalHydraulics-onePhase-blue

.. |tutorials/reactorCases/3D_genericParticleBedNTP-ThermalMechanics-extendedThermoMechanics-badge| image:: https://img.shields.io/badge/ThermalMechanics-extendedThermoMechanics-blue

.. |tutorials/reactorCases/3D_genericParticleBedNTP-Multiphysics-looseCoupling-badge| image:: https://img.shields.io/badge/Multiphysics-looseCoupling-orange

.. |tutorials/reactorCases/Godiva_SN-Neutronics-SN-badge| image:: https://img.shields.io/badge/Neutronics-SN-blue