.. _couplingGF:


Achieving coupled solutions 
===========================

Coupling logic
--------------

The coupling between physics is achieved by projecting coupling variables from
the mesh they are calculated, to the mesh they need to be used. The details of
the coupling can be specified in *system/regionsDict*. In the
sub-dictionary ``mappings``, for each region one can select the fields to map
*onto* it. This is done by creating a *subDict* named after the region *from*
which the fields are mapped. For instance, if a field needs to mapped into the
*fluidRegion* from the *neutroRegion*, the specifics of the mapping are found under
*regionsDict/mappings/fluidRegion/neutroRegion*. In this *subDict*, one can
specify the name of the field of the original mesh in the *sourceFields* entry
(e.g., ``powerDensity`` in the neutroRegion) and the name of the field onto which
the original field is mapped in the *targetFields* entry (e.g.
``powerDensityStructure`` in the fluidRegion). This routine is templated, hence
the user doesn't need to specify the field type (i.e. scalar or vector). A
detailed usage of this new coupling routine can be found in any multi-physics
tutorial, such as
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/>`_.

Example of loose coupling:

.. code :: cpp

    // In constant/regionsDict

    regionSolvers
    {
        Level_0
        {
            fluidRegion     onePhase;
            neutroRegion    diffusionNeutronics;
        }
    }

    mappings
    {
        fluidRegion
        {
            neutroRegion // to fluidRegion
            {
                sourceFields    ( powerDensity );
                targetFields    ( powerDensityNeutronics );
            }
        }
        neutroRegion
        {
            fluidRegion // to neutroRegion
            {
                sourceFields    ( T      thermo:rho );
                targetFields    ( TCool  rhoCool );
            }
        }
    }

Heterogeneous (boundary) vs Homogeneous (domain overlap) coupling
---------------------------------------------------------------

GeN-Foam can operate in two different
modes, depending on how the temperature of the structures is calculated and as
shown in the figure below:

- Homogeneous / domain overlap: the thermal-hydraulics solver is responsible for calculating temperatures throughout the physical domain. This is the most typical case, where structures are assumed to be treated as sub-scale structures in a porous-medium treatment.
- Heterogeneous / boundary coupling: the thermal-hydraulics solver and the thermal-mechanics solver are responsible for calculating temperatures in different parts of the domain. This could be used for instance when simulating a core with a porous-medium approach, and a large solid reflector using the thermo-mechanics solver. Or it could be used to simulate a fuel pin, assembly or entire code in a heterogeneous manner.

Of course, it is also possible to have hybrid approaches, where the temperature
in one structure is calculated in part based on the temperatures predicted in
the sub-scale structure by the thermal-hydraulics sub-solver, and partly by the
thermo-mechanical solver itself (where there is not overlap with the
thermo-hydraulics domain).

Tutorial `2D_fullCoupling <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/2D_fullCoupling>`_
has been created to allow users to play around with the couplings and understand
their logic.

.. figure:: ../../images/HetHom.png
    :width: 500
    :alt: Heterogeneous / Homogeneous meshes overlaps

    Heterogeneous / Homogeneous meshes overlaps.

.. note ::
    GeN-Foam can
    be used to map *any* scalar or vectorial field to *any* scalar of vectorial
    fields on a different mesh, enabling the simulation of any arbitrarily coupled multi-physics
    simulation that leverages the currently existing libraries.

.. note ::

    There is no need for the thermal-hydraulics, thermal-mechanics, and
    neutronics domains to be the same. GeN-Foam will project fields in a "clever" way
    whenever there is no overlap between 2 domains.


.. note ::

    It is possible not to solve for displacements in the thermo-mechanics
    sub-solver by setting to *false* the keyword *solveDisplacement* in
    *thermoMechanicalProperties*.

.. Note that temperature handling depends on mesh overlap:

.. - **Overlap with thermal-hydraulics mesh:**  
..   The solid temperature is forced to match the sub-scale structure temperature
..   predicted by thermal-hydraulics.

.. - **No overlap:**  
..   Temperature is obtained by solving a heat diffusion equation, using material
..   properties specified in ``thermoMechanicalProperties``.

.. - **Overlap with the neutronics mesh:**  
..   The region receives a volumetric heat source from the neutronics-calculated
..   power field.

Here is a list of the fields which are most commonly coupled across physics and
their names

.. list-table:: Common field correspondence table
    :widths: 50 50 50 50
    :header-rows: 1

    * - Physics Field
      - Thermal-hydraulics
      - Neutronics
      - Thermo-mechanics
    * - Coolant temperature
      - T
      - TCool
      - N/A
    * - Coolant velocity
      - U
      - U
      - N/A
    * - Coolant density
      - thermo:rho
      - rhoCool
      - N/A
    * - Fuel temperature (avg)
      - T.fuelAvForNeutronics
      - TFuel
      - TFuel
    * - Clad temperature (avg)
      - T.cladAvForNeutronics
      - TClad
      - TFuel (if linked fuel)
    * - Structures temperature
      - T.passiveStructure
      - TStruct
      - TStructFromTH/T\*
    * - Power density (structures)
      - powerDensityNeutronics
      - powerDensity
      - powerDensityNeutronics
    * - Power density (liquid)
      - powerDensityNeutronicsToLiquid
      - secondaryPowerDensity
      - N/A
    * - Displacement
      - N/A
      - disp
      - meshDisp/disp \*\*
    * - Coolant viscosity
      - mu
      - mu
      - N/A
    * - Coolant phase
      - alpha+phaseName e.g. alpha.liquid (twoPhase) alpha (onePhase)
      - alpha ***
      - N/A
    * - Thermal turbulent diffusivity
      - alphat (if exists)
      - alphat
      - N/A


\* In thermo-mechanics, the heat diffusion equation is solved only in the non-porous zones, while temperature from the thermal-hydraulics is expected in
the porous structures. Hence, the structures temperatures need to be mapped onto
the ``TStructFromTH`` field.

** The thermo-mechanics solver computes the displacement using the linear elastic
formulation. This result in the field named ``disp``. In order to compute the
neutronics deformation, the radial component of disp is combined to an axial
component computed as the thermal deformation of fuel/control rods (if
available). Hence, this is the field which is normally expected to be coupled to
the neutronics field named disp.

*** Most of the MSR cases which involve liquid fuel (and as a result the need to
map the coolant phase from TH to neutronics) are single-phase, hence simply
mapping alpha (TH) to alpha (neutronics) will work. However, in case a
multi-physics liquid-fuel, the coolant fraction that needs to be mapped is the
sum of ``alpha.liquid`` and ``alpha.vapour``. Such a field is currently not
existing in the thermal-hydraulics solver. However, using the conventional
OpenFOAM functions one can create a new field runTime and map it if needed. An
example of such approach is shown below.

In the *controlDict* add:

.. code :: cpp

    // In system/controlDict

    functions
    {
        newFieldCreation
        {
            type			coded;
            libs			("libutilityFunctionObjects.so");
            name			newFieldCreation;
            executeControl	timeStep;
            executeInterval 1;
            enabled			true;
            region          fluidRegion; //or whatever other region, but important to specify otherwise it is defaulted to region0 and the fields are not found

            codeExecute
            #{
                // Lookup the fields needed for the creation of the new field
                const volScalarField& field1 = mesh().lookupObjectRef<volScalarField>("field1"); // Use correct names of the fields
                const volScalarField& field2 = mesh().lookupObjectRef<volScalarField>("field2");

                // Creation of the new field

                static autoPtr<volScalarField> newField;

                if (!mesh().foundObject<volScalarField>("newField"))
                {
                    newField.reset
                    (
                        new volScalarField
                        (
                            IOobject
                            (
                                "newField",
                                mesh().time().timeName(),
                                mesh(),
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                            ),
                            field1 + field2 // here as an example the sum of the two fields is considered
                        )
                    );
                }
                else
                {
                    // Insert here any operation required
                    newField() = field1+field2;
                    newField().correctBoundaryConditions();
                }
            #};
        }
    }


The *loops*
-----------

Many multi-physics simulations might require large flexibility on the
time-loops. For instance, while some physics might be tightly coupled, others
are only loosely coupled. This is, for instance, often the case for multi-scale
simulations, where physics are tightly coupled within a scale and loosely across
scales. In order to create this additional flexibility, dedicated "loops" are created. For example, the ``picardLoop`` consists of
a list of tightly coupled physics. During run-time, when the physics are
corrected, if they are part of a ``picardLoop``, the physics will be
corrected iteratively, until a user-input convergence criterion is met.
*picardLoops* are also defined in *system/regionsDict* in
the ``regionSolvers`` sub-dictionary. Each entry of the subDict corresponds
to a different ``picardLoop``, characterized by a list of subSolvers (i.e.
the solvers that are tightly coupled), a minimum residual to reach before the
Picard iterations are interrupted and a maximum number of iterations. Note that
*picardLoops* can be nested within one other, possibly leading to
tree-like structure as shown below.

.. mermaid::

    graph TD
        A[GeN-Foam] --> B[Solver 1]
        A --> C[MultiPhysicsSolver 1]
        C --> D[Solver 2]
        C --> E[MultiPhysicsSolver 2]
        E --> F[Solver 3]
        E --> G[Solver 4]


Which can be translated in *system/regionsDict*:

.. code :: cpp

    // In constant/regionsDict

    regionSolvers
    {
        Level_0
        {
            Solver1                 physicsB;
            MultiPhysicsSolver1     picardLoop;
        }
        MultiPhysicsSolver1
        {
            solvers
            {
                Solver2                 physicsD;
                MultiPhysicsSolver2     picardLoop;
            }
            minResidual     0.00005;
            maxIterations   3;
        }
        MultiPhysicsSolver2
        {
            solvers
            {
                Solver3     physicsF;
                Solver4     physicsG;
            }
            minResidual     0.00005;
            maxIterations   3;
        }
    }


Note that all the physics that are specified in the ``regionSolvers`` dict are
loosely coupled, whereas only those specified within a ``picardLoop`` are
tightly coupled.

The previous example would correspond to a time loop as shown below.

.. mermaid::

    graph TD
        S@{ shape: sm-circ, label: "Small start" } --> A
        A[New time step] --> B[Solver 1]
        B --> C[Solver 2]
        C --> D[Solver 3]
        D --> E[Solver 4]
        E --> F{Is MultiPhysicsSolver 2 converged?}
        F -- No --> D
        F -- Yes --> G{Is MultiPhysicsSolver 1 converged?}
        G -- No --> C
        G -- Yes --> A
        G --> End@{ shape: framed-circle, label: "Stop" }

Current loops available in *GeN-Foam* are:

toctreeHere




Typical coupling logic for a nuclear reactor
--------------------------------------------

The standard sequence used to solve multi-physics nuclear reactor problems is
summarized with reference to the figure below.

.. figure:: ../../images/couplingFields.png
    :width: 500
    :alt: Typical coupling for logic porous-medium  reactor calculations

    Typical coupling for logic porous-medium  reactor calculations.

Given an input volumetric power density :math:`q`, the thermal-hydraulics
sub-solver computes the resulting fields:

- fluid temperature :math:`T`,
- fluid density :math:`\rho`,
- fluid velocity :math:`U`,
- and structural temperatures :math:`T_s`.

For two-phase treatments, the fields :math:`T`, :math:`\rho`, and :math:`U`
represent mass-weighted mixture values. The velocity field :math:`U` is only
used for solver coupling in molten salt reactors (MSRs), where it is required to
advect delayed neutron precursors.

The power density :math:`q` may refer either to a sub‑scale structure
(e.g. a fuel rod) or to the fluid itself in liquid-fuel systems, depending on the
reactor type. The symbol :math:`T_s` collectively denotes temperatures of solid
structures—fuel, cladding, control rod drivelines, wrapper tubes, the diagrid,
and others—depending on which structural models are assigned to each cell zone.

The neutronics sub-solver computes the updated fuel volumetric power density
:math:`q` based on the coupling fields. Depending on the neutronics model
(diffusion, :math:`S_N`, or :math:`SP_3`), the solver may incorporate feedbacks
from:

- coolant temperature :math:`T` and density :math:`\rho`,
- average fuel and cladding temperatures :math:`T_s`,
- fuel axial expansion and core radial expansion, collectively :math:`d`,
- temperatures predicted by the thermal-mechanics sub-solver.

Whenever macroscopic cross-sections are parameterized against these fields, the
solver can account for the associated reactivity feedbacks. For point-kinetics,
feedbacks follow standard coefficient-based models.

The thermal-mechanics sub-solver computes temperatures in solid regions and an
overall displacement field used to deform the neutronics mesh. The displacement
is decomposed into axial and radial components, collectively denoted
:math:`d`, and passed to the neutronics module to model expansion-related
feedback.


