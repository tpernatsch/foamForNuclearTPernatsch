.. _userguide_neutronics_initialAndBC:


Initial and boundary conditions
-------------------------------

Initial and boundary conditions (IC and BC) adopt the usual OpenFOAM logic. 
OpenFOAM provides most of the boundary conditions one may need for
neutronics models. In addition, a few boundary conditions have been
included in the neutronics module of foamForNuclear:

.. toctree::
   :maxdepth: 1
   :glob:

   ../../../cppapi/generated/fvPatchFields/neutronics/**

.. raw:: html

   <br><br>



For the initial conditions of fluxes, 
the user can either specify different IC and BC for each one of the energy
groups (with fluxes that must be named *fluxStar0*, *fluxStar1*, etc…) or
provide the same IC and BC to all fluxes by using the *defaultFlux* field. In
the case of SP3 calculations, the IC and BC for the second moment can be imposed
either for each energy (using fields named *fluxStar20*, *fluxStar21*, etc…),
or to all energies by using the *defaultFlux2* field. When both *defaultFlux*
and *fluxStar...* are present, the solver gives priority to *fluxStar<g>*. In
the case of SN calculations, it is suggested not to modify the boundary
conditions and to use the *defaultFlux* file (an example is provided in the
`Godiva_SN <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/Godiva_SN/0/neutroRegion/defaultFlux>`_ tutorial). 
When employing the adjoint solver, the user will have to add the
fields *adjointDefaultPrec* and *adjointDefaultFlux* in the initial time.


IC and BC for precursors do not have to be specified for standard reactors. On
the other hand, they should be specified in the case of liquid fuel reactors
(e.g., Molten Salt Reactors). This is possible by creating a *defaultPrec*
field, in case the same conditions apply to all precursor groups, or by creating
the fields named *prec0*, *prec1*, etc., in case different conditions must be
provided for different precursor groups.

.. note ::

    Boundary conditions must be applied to ``fluxStar<g>`` and not to ``flux<g>``
    since the neutronics module solves for these variables. ``fluxStar<g>`` represent continuous
    fluxes, while ``flux<g>`` represent the real fluxes. They differ only in case
    discontinuity factors are employed (see :ref:`FIORINA2016212 <FIORINA2016212>`).





