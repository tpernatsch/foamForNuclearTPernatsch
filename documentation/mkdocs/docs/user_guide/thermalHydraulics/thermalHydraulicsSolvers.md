# Sub-solvers

Thermal-hydraulics calculations are performed by classes derived from *thermalHydraulicsModel.H* that contain specific sub-solvers:

- *onePhase* for single-phase calculations, using the formulation proposed in Refs. [@Radman2019ADesign] [@RADMAN2021111178] [@RADMAN2021111422] (see [*onePhase.H*](../../classes/thermalHydraulics/solvers/onePhase/onePhase.md))
- *twoPhase* for adjoint diffusion calculations, using the formulation proposed in Refs. [@Radman2019ADesign] [@RADMAN2021111178] [@RADMAN2021111422] (see [*twoPhase.H*](../../classes/thermalHydraulics/solvers/twoPhase/twoPhase.md))

For the user, the derived classes translate into runtime selectable models. The specific sub-solver to be used in a simulation can be selected in the solvers dictionary like explained in the coupling section.
