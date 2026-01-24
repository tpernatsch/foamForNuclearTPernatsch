.. _userguide_neutronics_settingcase:

Setting a case
--------------

Initial and boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As in all standard OpenFOAM solvers, initial conditions (IC) and boundary conditions
(BC) should be provided in the "0" folder, or the folder corresponding to the
``startTime`` of the simulation, if different than 0. In the case of neutronics,
the user can either specify different IC and BC for each one of the energy
groups (with fluxes that must be named *fluxStar0*, *fluxStar1*, etc…) or
provide the same IC and BC to all fluxes by using the *defaultFlux* field. In
the case of SP3 calculations, the IC and BC for the second moment can be imposed
either for each energy (using fields named *fluxStar20*, *fluxStar21*, etc…),
or to all energies by using the *defaultFlux2* field. When both *defaultFlux*
and *fluxStar...* are present, the solver gives priority to *fluxStar...*. In
the case of SN calculations, it is suggested not to modify the boundary
conditions and to use the *defaultFlux* file (an example is provided in the
Godiva_SN tutorial). When employing the adjoint solver, you will have to add the
fields *adjointDefaultPrec* and *adjointDefaultFlux* in your initial time.

In addition to the standard OpenFOAM BC, an albedo boundary condition (see
:ref:`albedoSP3FvPatchField.H <albedoSP3FvPatchField>`) is available in GeN-Foam
for diffusion and SP3 calculations and can be used according to the following
syntax:

.. code :: cpp

    type            albedoSP3;
    gamma           0.5;            // defined as (1-alpha)/(1+alpha)/2, alpha being the albedo coefficient
    diffCoeffName   Dalbedo;        // not to be changed
    fluxStarAlbedo  fluxStarAlbedo; // not to be changed
    forSecondMoment false;          // true in case it is a condition for a second moment flux (for SP3 calculations)
    value           uniform 1;


Please note that the boundary condition needs to be set both for the first and
second moments in SP3.

IC and BC for precursors do not have to be specified for standard reactors. On
the other hand, they should be specified in the case of liquid fuel reactors
(e.g., Molten Salt Reactors). This is possible by creating a *defaultPrec*
field, in case the same conditions apply to all precursor groups, or by creating
the fields named *prec0*, *prec1*, etc., in case different conditions must be
provided for different precursor groups.

.. note ::

    Boundary conditions must be applied to ``fluxStar...`` and not to ``flux...``
    since GeN-Foam solves for these variables. ``fluxStar...`` represent continuous
    fluxes, while ``flux...`` represent the real fluxes. They differ only in case
    discontinuity factors are employed (see :ref:`FIORINA2016212 <FIORINA2016212>`).


Setting the weighting in point-kinetics calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A correct evaluation of the reactivity worth of delayed neutron precursors in
MSRs, as well as of the impact of temperatures on reactivities, normally
requires the knowledge of the adjoint flux. In GeN-Foam, the field
*oneGroupFlux* is used by the point kinetic solver for weighting temperatures,
densities and precursors. When fluxes are not calculated via a spatial
neutronics calculation, one has to manually provide the *oneGroupFlux* in
*0/neutroRegion*. As an alternative, one can use the *initialOneGroupFluxByZone*
keyword in *nuclearData* (see `1D_MSR_pointKinetics
<https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/1D_MSR_pointKinetics/rootCase/constant/neutroRegion/nuclearData>`_). Please notice that:

- If calculated fluxes are available in *neutroRegion*, these will be user to recalculate and overwrite *oneGroupFlux*.
- If no fluxes are available, the neutronics sub-solver will use the provided *oneGroupFlux*
- If the *initialOneGroupFluxByZone* keyword is used in *nuclearData*, this will be used to overwrite *oneGroupFlux*


Discretization and solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Details for discretization and solution of equations are handled in a standard
OpenFOAM way, i.e., through the *fvSolution* and *fvSchemes* dictionaries in
*system/neutroRegion*.
