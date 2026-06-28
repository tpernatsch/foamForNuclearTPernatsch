Fmu Cases Tutorials
===================

This section provides tutorial-style case studies showing how foamForNuclear can be coupled to an external FMU-based balance of plant and how control logic can be incorporated into the coupled simulation workflow. The first subsection demonstrates a 2D-wedge ALFRED model with an FMU representing the balance of plant, establishing the baseline coupling between reactor physics and external system dynamics. The second subsection extends this idea by coupling a 2D point-kinetics model to an FMU, illustrating a reduced-order reactor model within the same FMU interaction pattern. The final subsection then focuses on a power–temperature–momentum controller, showing how feedback control can be layered on top of the coupled reactor–FMU setup to regulate system behavior.

.. toctree::
    :maxdepth: 1

    2D_LFRpowerPlant
    2D_PKCoupleFMI
    powerTemperatureMomentumControl
