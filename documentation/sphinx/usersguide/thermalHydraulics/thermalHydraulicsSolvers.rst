===========
Sub-solvers
===========

Thermal-hydraulics calculations are performed by classes derived from
:ref:`thermalHydraulicsModel.H <thermalHydraulicsModel>` that contain specific
sub-solvers:

:onePhase: for single-phase calculations, using the formulation proposed in
    Refs. [@Radman2019ADesign] [:ref:`RADMAN2021111178 <RADMAN2021111178>`]
    [:ref:`RADMAN2021111422 <RADMAN2021111422>`] (see :ref:`onePhase.H <onePhase>`)
:twoPhase: for adjoint diffusion calculations, using the formulation proposed in
    Refs. [@Radman2019ADesign] [:ref:`RADMAN2021111178 <RADMAN2021111178>`]
    [:ref:`RADMAN2021111422 <RADMAN2021111422>`] (see :ref:`twoPhase.H <twoPhase>`)

For the user, the derived classes translate into runtime selectable models. The
specific sub-solver to be used in a simulation can be selected in the solvers
dictionary like explained in the coupling section.
