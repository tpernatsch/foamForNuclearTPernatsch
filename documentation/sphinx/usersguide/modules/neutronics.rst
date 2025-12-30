.. _usersguide_neutronics:


Neutronics
==========


Introduction
------------

The neutronics class (see :ref:`neutronics.H <neutronics>`) is a high-level
class that contains essential data and variables that are common to various
neutronics models. In articular, the neutronics class handles the variables that
are included in the *0/uniform/reactorState*,
*constant/neutroRegion/nuclearData* and
*constant/neutroRegion/neutronicsProperties* dictionaries.


.. _neutronics_the-reactorstate-dictionary:

The *reactorState* dictionary
-----------------------------

The *reactorState* dictionary is found under the *timeFolder/uniform/*
sub-folder. It includes essentially 3 keywords:

:keff: is used in the spatial kinetics solvers as an initial guess for keff when
       doing an eigenvalue calculation. It is then updated automatically at each
       time step (i.e., at each power iteration) with the calculated value of
       keff. When performing a transient calculation with the spatial kinetics
       solvers, *keff* is instead used to divide the neutron source term and is
       not updated during the simulation. Typically, to run spatial kinetics
       transient simulations, one first runs an eigenvalue calculation. The
       resulting *keff* will be the one that makes the reactor critical in a
       subsequent transient simulation. *keff* is disregarded by the point
       kinetics sub-solver.
:pTarget: is used in the spatial kinetics solvers as target power when doing an
          eigenvalue calculation. It is also used by the point kinetics
          sub-solver, but only to correctly plot results. As power, GeN-Foam
          uses what it finds under powerDensity of the neutroRegion, or under
          the powerDensity of the fluidRegion if it does not find a
          powerDensity in the neutroRegion. To correctly plot point kinetics
          results, pTarget must be consistent with the mentioned power
          densities.
:precursorPowers: can be read by the point kinetics sub-solver in case the user
                  wishes to set initial concentrations of precursors. If not
                  found, precursor concentrations are initialized to be in
                  equilibrium with the starting conditions (i.e. a steady state
                  is assumed).


All GeN-Foam neutronics models can be used for liquid-fuel reactors. One can
activate this option using the *liquidFuel* keyword in */system/controlDict*.
Of course, in such cases, one should pay attention to setting proper boundary
conditions for the precursors.

A commented *reactorState* can be found in `3D_SmallESFR
<https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/0/uniform/reactorState>`_.

NB: Please note that in parallel calculations, the updated *reactorState* can be found in *timeStep/uniform/*.


IMPORTANT: The powerDensity file written to disk is the power density calculated
by neutronics (sigmaPowers multiplied by fluxes), DIVIDED by the fuelFractions
indicated in nuclearData. This means that it provides the power density IN the
fuel, not spread over the cross-section homogenization region. For instance, if
you have an assembly with its own one-group cross-section set and you specify
that the fuel fraction is 0.3, powerDensity will be equal to:

.. math::
    q''' = \frac{\kappa \Sigma_f \phi}{\alpha_{fuel}} = \frac{\text{sigmaPower} \times \text{flux}}{\text{fuelFraction}}



Models
------

Neutronics calculations are performed by classes derived from *neutronics* that
contain specific sub-solvers:

:pointKinetics: Point-kinetics (:ref:`pointKineticNeutronics.H <pointKineticNeutronics>`)
:diffusionNeutronics: Diffusion (:ref:`diffusionNeutronics.H <diffusionNeutronics>`)
:adjointDiffusion: Adjoint diffusion (:ref:`adjointDiffusionNeutronics.H <adjointDiffusionNeutronics>`)
:SP3Neutronics: Diffusion in :math:`SP_3` (:ref:`SP3Neutronics.H <SP3Neutronics>`)
:SNNeutronics: Discrete ordinates (:math:`S_N`) (:ref:`SNNeutronics.H <SNNeutronics>`)

For the user, the derived classes translate into runtime selectable models. The
specific sub-solver to be used in a simulation can be selected at runtime in the
*constant/neutroRegion/neutronicsProperties* dictionary.

The choice of the model is achieved by selecting the wanted solver in 
*regionSolvers* depending on whether the neutronics
solvers need to be part of a tightly coupled loop or not (see
:ref:`Coupling solvers <userguide_coupling_the-controlDict-dictionary>`).

*adjointDiffusion* has been developed only as an eigenvalue solver. The others
can be used for transient calculations. However, the SN transient solver has not
been tested. In addition, it is currently not accelerated, thus extremely slow
(it can require hundreds of iterations per time step).


The *neutronicsProperties* dictionary
-------------------------------------

The *neutronicsProperties* dictionary is found under *constant/neutroRegion/*
and it can be used to set the type of neutronics simulation by using the
following keywords:

:eigenvalueNeutronics: should be set to ``true`` for eigenvalue calculations,
                       false for transients.
:externalSourceNeutronics: should be set to ``true`` for external neutron source
                           calculations. The ``eigenvalueNeutronics`` variable
                           should be put to ``false`` and the ``keff = 1``.

One can find detailed, commented examples in most tutorials. See for instance
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/neutronicsProperties>`_ (single phase).

N.B: The parameter *model* used to define what type of simulation needs to be
performed as been replaced by the selection of model in the *system/regionsDict*
(see :ref:`Coupling solvers <userguide_coupling>`).


.. _userguide_neutronics_nuclearData:


The *nuclearData* dictionary
----------------------------

In GeN-Foam, cross-sections and several other neutronics properties are handled
by the :ref:`XS.H <XS>` class. Detailed explanations on the file format are provided in
:ref:`XS.H <XS>` and in the tutorials (e.g `3D_SmallESFR
<https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/nuclearData>`_).

The *nuclearData* dictionary can be found under *constant/neutroRegion/*. It
contains all basic nuclear properties for the reference and perturbed reactor
states. For instance, including ``TFuel`` in the ``reference`` state and a perturbed state
represents the temperatures at which the reference and perturbed cross-sections
have been calculated, respectively. Radial Basis Function interpolation is
performed by GeN-Foam between reference and perturbed reactor states. It is
possible to provide the XS set with multiple perturbation variables (see example
below). If no perturbed state data are provided, the reference cross-sections
are used.

Special field for axial and radial expansions are provided as ``axExp`` and
``radExp`` (see `3D_SmallESFR
<https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/nuclearData>`_).

Nuclear data can be generated using any nuclear code.

:`serpentToFoam <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tools/serpentToFoam/serpent2.1.23>`_:
    routines provided with GeN-Foam (in the *tools* folder) is an Octave script
    that automatically converts Serpent output files into the nuclear data files
    employed by GeN-Foam.
:`openmcToFoam <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tools/openmcToFoam>`_:
    Python package provided with GeN-Foam automatically converts OpenMC output
    into nuclear data files.

The entry ``discFactor`` is used only if discontinuity factors have to be used.
The term ``integralFlux``, is used only if the automatic adjustment of
discontinuity factors is performed (see :ref:`FIORINA2016212 <FIORINA2016212>`). Nonetheless, these entries
should always be present.



XS parametrization
------------------

GeN-Foam features XS parametrization using the Radial Basis Function
interpolation scheme on any field provided by GeN-Foam. This method allows to
interpolate the XS using multiple parameters/perturbations (see example below).

It is possible to select different radial basis function based on the
polyharmonic splines using the ``polyharmonicSplineMode`` keyword. The figure
below shows the radial basis function influence on the interpolation. The
default mode is ``1``, which guarantee a linear interpolation:

:`1`: :math:`\phi(r) = |r|`
:`2`: :math:`\phi(r) = r^2 \ln(r)`
:`3`: :math:`\phi(r) = |r^3|`
:`4`: :math:`\phi(r) = r^4 \ln(r)`

.. figure:: ../images/fig_RBF_interpolationOrders.png
    :width: 500
    :alt: Effects of Polyharmonic Spline Radial Basis Function order on arbitrary set of XS points.

    Effects of Polyharmonic Spline Radial Basis Function order on arbitrary set of XS points.

In combination to the RBF interpolation, three laws are currently provided to
the user to modify the behavior and improve interpolation accuracy:

:``lin``: linear
:``sqrt``: square root
:``log``: logarithmic

It is possible to assign the law through the ``xsVariables`` sub dictionary in
*nuclearData* with the name of the field. If one or several fields provided in
``xsVariables`` are not default to GeN-Foam (e.g ``Tmatrix``), the code will
automatically create it in the neutronics region and can be used for additional
coupling with other solvers (see the :ref:`coupling page <userguide_coupling>`).

.. code :: cpp

    xsVariables
    {
        TFuel       log;
        rhoCool     lin;
    }

    states
    (
        reference // Mandatory name not to be modified
        {
            TFuel   900;
            rhoCool 4125;

            zones
            (
                zone1
                {
                    ...
                }
                ...
            );
        }

        Tfuel1200K // Arbitrary name
        {
            TFuel   1200;
            #include "XSTfuel1200K" // OpenFOAM shortcut to attach file content at this location
        }

        rhoCool3500kgm3
        {
            rhoCool 3500;
            #include "XSrhoCool3500kgm3"
        }

        Tfuel1200KandRhoCool3500kgm3
        {
            TFuel   1200;
            rhoCool 3500;
            #include "XSTfuel1200KandRhoCool3500kgm3"
        }
    );


.. note ::

    Cross-sections must be expressed according to the International System of
    Units (so m, not cm).

.. note ::

    ``defaultPrec`` has 1/m3 units except for the adjoint solver that needs 1/m2/s.

.. note ::

    The *nuclearData* file must always be present, even when not parametrizing
    cross-sections. If no parametrization is needed, the ``zones`` card must be
    left "blank" as:


.. code :: cpp

    xsVariables
    {}

    states
    (
        reference
        {
            zones
            ();
        }
    );


One can find more details on all the parameters in the :ref:`XS.H <XS>` file and commented
examples of *nuclearData* in the tutorials
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/nuclearData>`_ (for diffusion or SP3),
`Godiva_SN <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/Godiva_SN/constant/neutroRegion/nuclearData>`_ (for discrete ordinates) and
`2D_onePhaseAndPointKineticsCoupling <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/2D_onePhaseAndPointKineticsCoupling/rootCase/constant/neutroRegion/nuclearData>`_ (for point kinetics).
`2D_onePhaseAndSubcriticalPointKineticsCoupling <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/2D_onePhaseAndSubcriticalPointKineticsCoupling/rootCase/constant/neutroRegion/externalSource>`_ (for subcritical point kinetics).




.. _userguide_neutronics_quadratureSet:

The *quadratureSet* dictionary
------------------------------

An additional dictionary is needed to provide the quadrature set when performing
discrete ordinate calculations. The *quadratureSet* dictionary is found under
*constant/neutroRegion/*. It contains the quadrature set for discrete ordinate
calculations. One can find examples of three different quadrature sets in the
tutorial `Godiva_SN
<https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/Godiva_SN/constant/neutroRegion/>`_.
S4 and S8 Chebyshev-Legendre quadrature sets can be found in `Godiva_SN
<https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tools/chebichevLegendreQuadratureSets/>`_.


.. _userguide_neutronics_crmove:

The *CRMove* dictionary
-----------------------

The *CRMove* dictionary can be found under *constant/neutroRegion/*. It contains
input data for control rod movement. Control rods can be moved from the initial
position to a new one by selecting the initial and final time of the
insertion/extraction and the speed of insertion/extraction (positive speed for
insertion).

One can find a commented example in the tutorial
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/CRmove>`_, though this option is not actually used in the tutorial.


Initial and boundary conditions
-------------------------------

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
----------------------------------------------------

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
---------------------------

Details for discretization and solution of equations are handled in a standard
OpenFOAM way, i.e., through the *fvSolution* and *fvSchemes* dictionaries in
*system/neutroRegion*.


--------------------------
Subcritical point-kinetics
--------------------------

To use the subcritical point-kinetics, the user has to add the
*constant/neutroRegion/externalSource*. The file contains a flag to activate the
external neutron source (``isExternalSource``).

Several parameters related to a spallation source are included such as the
energy per source particle in J/source particle and the neutron yield of the
reaction in neutrons/source particle*

An external source modulation timetable is provided to manually modulate the
source strength.

In the case of an FMI coupling, it is possible to use the
``externalSourceModulationNameFromFMU`` entry to change the external source
modulation through an FMI. To use it, the mode must be ``transient``.


.. note ::
    ``defaultPrec`` field has 1/m3 units except for the adjoint solver that needs 1/m2/s.

.. note ::
    In point kinetics, pTarget (in reactorState) MUST be the same as the one used for reaching the steady-state. As power, GeN-Foam uses what it finds under powerDensity, or under the powerDensity of the fluidRegion if it does not find a powerDensity in the neutroRegion. pTarget does not enter the calculation, it is used simply to plot the results