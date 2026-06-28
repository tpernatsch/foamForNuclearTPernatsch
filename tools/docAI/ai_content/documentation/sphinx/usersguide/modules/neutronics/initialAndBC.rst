.. _userguide_neutronics_initialAndBC:


Initial and boundary conditions
-------------------------------

Initial and boundary conditions (IC and BC) follow the usual OpenFOAM logic. OpenFOAM provides most of the boundary conditions needed for neutronics models. In addition, a few boundary conditions have been included in the neutronics module of foamForNuclear:

.. toctree::
   :maxdepth: 1
   :glob:

   ../../../cppapi/generated/fvPatchFields/neutronics/**

.. raw:: html

   <br><br>



For the initial conditions of fluxes, the user can either specify different IC and BC for each energy group (with fluxes named *fluxStar0*, *fluxStar1*, etc.) or provide the same IC and BC to all fluxes by using the *defaultFlux* field. For SP3 calculations, the IC and BC for the second moment can be imposed either for each energy (using fields named *fluxStar20*, *fluxStar21*, etc.) or for all energies by using the *defaultFlux2* field. When both *defaultFlux* and *fluxStar...* are present, the solver gives priority to *fluxStar<g>*.

For SN calculations, it is suggested not to modify the boundary conditions and to use the *defaultFlux* file (an example is provided in the `Godiva_SN <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/Godiva_SN/0/neutroRegion/defaultFlux>`_ tutorial). When employing the adjoint solver, the user must add the fields *adjointDefaultPrec* and *adjointDefaultFlux* at the initial time.

IC and BC for precursors do not have to be specified for standard reactors. On the other hand, they should be specified for liquid fuel reactors (e.g., Molten Salt Reactors). This can be done by creating a *defaultPrec* field, if the same conditions apply to all precursor groups, or by creating the fields named *prec0*, *prec1*, etc., if different conditions must be provided for different precursor groups.

.. note ::

    Boundary conditions must be applied to ``fluxStar<g>`` and not to ``flux<g>``
    since the neutronics module solves for these variables. ``fluxStar<g>`` represent continuous
    fluxes, while ``flux<g>`` represent the real fluxes. They differ only when discontinuity factors are employed (see :ref:`FIORINA2016212 <FIORINA2016212>`).