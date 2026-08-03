.. _convertingOpenFoamSolver:

Converting a native OpenFOAM solver into a foamForNuclear class
=================================================================

This tutorial explains, step by step, how to transpose a standard
(``main()``-based) OpenFOAM solver application into a run-time selectable
foamForNuclear/GeN-Foam solver *class*. It is written from a direct comparison
of:

- ``applications/modules/openFoamImported/compressibleInterFoamNATIVE``: an
  unmodified copy of the OpenFOAM ``compressibleInterFoam`` application, kept
  in the repository purely as a reference for what the "before" state looks
  like (it is *not* compiled — there is no top-level ``Make`` entry for it),
  and
- ``applications/modules/openFoamImported/compressibleInterFoam``: the same
  solver, transposed into the foamForNuclear class format and actually
  compiled into ``libOFSolvers``.

A second, independent example of the same conversion (OpenFOAM's
``rhoPimpleFoam``) can be found in
``applications/modules/thermalHydraulics/onePhaseRhoPimpleFoam``, and is worth
reading side-by-side with this tutorial since it applies the exact same
pattern with fewer leftover inconsistencies.

.. note ::
    The disclaimer carried in every converted file in this repository applies
    here too: *"The source code may not be in the OpenFOAM coding style, and it
    might not be making use of inheritance of classes to full extent."* The
    goal of the conversion is to make the solver runtime-selectable and
    embeddable in a multi-region simulation, not to rewrite it as idiomatic,
    fully encapsulated C++.


Why this is needed
-------------------

A standard OpenFOAM solver is a ``main()`` function: it owns the whole
time loop, reads its fields as free-standing local variables (via literal
``#include`` of files such as ``createFields.H``), and is built as a
standalone executable (``EXE``). foamForNuclear/GeN-Foam instead runs an
arbitrary number of coupled regions in a single process
(see :ref:`Achieving coupled solutions <couplingGF>`), where each region is
solved by an object selected at run time from the *system/regionsDict*
dictionary and driven by a shared, outer time loop
(``applications/solvers/GeN-Foam/GeN-Foam.C``):

.. code :: cpp

    regionSolvers solvers(runTime);
    ...
    while (runTime.run())
    {
        ...
        forAll(solvers, i)
        {
            solvers[i].deformMesh();
            solvers[i].correctPhysics();
            solvers[i].correctBaffleLessFields();
        }
        adjustDeltaT(runTime, solvers);
        ...
    }

For a solver to be one of these ``solvers[i]`` objects, it must be
transposed into a class deriving from the abstract base class
``Foam::solver`` (``src/multiRegion/solver/solver.H``). That is precisely what
turns ``compressibleInterFoamNATIVE`` into ``compressibleInterFoam``.


The base class: ``Foam::solver``
---------------------------------

Every wrapped solver derives, directly or indirectly, from
``src/multiRegion/solver/solver.H``:

.. code :: cpp

    class solver
    :
        public regIOobject
    {
    protected:
        dynamicFvMesh& mesh_;
        const Time&    runTime_;
        bool steady_;
        bool LTS_;

    public:
        TypeName("solver");

        declareRunTimeSelectionTable
        (
            autoPtr,
            solver,
            dynamicFvMesh,
            (dynamicFvMesh& mesh),
            (mesh)
        );

        solver(dynamicFvMesh& mesh);

        static autoPtr<solver> New(const word& solverName, dynamicFvMesh& mesh);

        virtual scalar maxDeltaT() = 0;
        virtual void   correctPhysics() = 0;
        virtual void   correctTightlyCoupledPhysics() { this->correctPhysics(); }
        virtual void   correctBaffleLessFields() {};
        virtual void   deformMesh() {};
        virtual scalar getResidual() = 0;

        inline virtual const Time&   runTime() { return runTime_; }
        inline virtual const fvMesh& mesh()    { return mesh_; }
    };

The important consequences for the conversion:

- The constructor takes a single ``dynamicFvMesh&`` — this replaces
  ``createMesh.H`` / ``createTime.H``: the mesh and time are handed to you.
- ``mesh_`` and ``runTime_`` already exist as protected members, so every
  bare ``mesh`` / ``runTime`` in the original solver becomes ``mesh_`` /
  ``runTime_`` (or the accessors ``mesh()`` / ``runTime()``).
- Three methods are pure virtual and **must** be implemented:
  ``correctPhysics()``, ``maxDeltaT()``, ``getResidual()``.
- The class is looked up at run time by name through a run-time selection
  table keyed on ``dynamicFvMesh&`` — this is what lets
  ``system/regionsDict`` refer to your solver by a plain word.


Step 1 — lay out the module directory
--------------------------------------

Keep a pristine, uncompiled copy of the native solver next to the converted
one — this is exactly what ``compressibleInterFoamNATIVE`` is for, and it is
invaluable for diffing later once you start modifying equation files.

Directory layout used by ``compressibleInterFoam``:

.. code :: text

    applications/modules/openFoamImported/
    ├── Make/
    │   ├── files          # shared by all solvers in this module directory
    │   └── options
    ├── compressibleInterFoamNATIVE/   # reference copy only, not built
    │   ├── compressibleInterFoam.C
    │   ├── createFields.H
    │   ├── UEqn.H / TEqn.H / pEqn.H / ...
    │   └── Make/
    └── compressibleInterFoam/          # the converted class
        ├── compressibleInterFoam.H
        ├── compressibleInterFoam.C
        ├── compressibleInterFoam.yaml  # doc metadata (see step 7)
        └── include/
            └── equations/
                ├── UEqn.H
                ├── TEqn.H
                ├── pEqn.H
                ├── alphaControls.H
                ├── alphaCourantNo.H
                ├── alphaEqn.H
                ├── compressibleAlphaEqnSubCycle.H
                ├── createAlphaFluxes.H
                ├── rhofs.H
                ├── alphaSuSp.H
                ├── setDeltaT.H
                └── setRDeltaT.H

Note that the converted version copies in *every* included ``.H`` file it
needs (even ones that in stock OpenFOAM live in a shared sibling directory
such as ``$FOAM_SOLVERS/multiphase/VoF``), so that the module is
self-contained and does not depend on files outside the repository.


Step 2 — write the class header
---------------------------------

Create ``<solverName>/<solverName>.H``. Its skeleton, taken from
``compressibleInterFoam.H``:

.. code :: cpp

    #ifndef compressibleInterFoam_H
    #define compressibleInterFoam_H

    #include "fvCFD.H"
    // ... every #include that the original solver's main()/createFields.H used
    #include "solver.H"

    namespace Foam
    {
    namespace solvers
    {

    class compressibleInterFoam
    :
        public solver
    {
    private:

        pimpleControl pimple_;

        volScalarField p_rgh_;
        volVectorField U_;
        surfaceScalarField phi_;
        twoPhaseMixtureThermo mixture_;
        volScalarField& alpha1_;
        // ... one member per field/object created by createFields.H

        fv::options& fvOptions_;
        IOMRFZoneList MRF_;

    public:

        TypeName("compressibleInterFoam");

        compressibleInterFoam(dynamicFvMesh& mesh_);
        virtual ~compressibleInterFoam() {}

        virtual void   correctPhysics();
        virtual void   correctTightlyCoupledPhysics();
        virtual scalar maxDeltaT();
        virtual scalar getResidual() { return scalar(0); }
        virtual void   correctBaffleLessFields() {}
        virtual void   deformMesh() {}

        // Optional finer-grained hooks (left empty if unused)
        void correctEnergy();
        void correctFluidMechanics();
        void correctCourant();
        void correctContErr();
        void printContErr();
        void calcCumulContErr();
    };

    } // End namespace solvers
    } // End namespace Foam

    #endif

Two rules drive this step:

1. **Namespace.** Wrapped solver classes live in ``Foam::solvers`` (plural),
   distinct from the base class ``Foam::solver`` (singular) which lives
   directly in ``Foam``.
2. **Member-variable conversion.** Every field, object or dictionary entry
   declared as a local/global variable in ``createFields.H`` (and any
   auxiliary include such as ``createAlphaFluxes.H``, ``createMRF.H``,
   ``createFvOptions.H``) becomes a **private member with a trailing
   underscore**. For example:

   .. list-table:: createFields.H → class members
       :widths: 50 50
       :header-rows: 1

       * - Native (``createFields.H``, global)
         - Wrapped (class member)
       * - ``volScalarField p_rgh(...)``
         - ``volScalarField p_rgh_;``
       * - ``volVectorField U(...)``
         - ``volVectorField U_;``
       * - ``twoPhaseMixtureThermo mixture(U, phi);``
         - ``twoPhaseMixtureThermo mixture_;``
       * - ``compressibleInterPhaseTransportModel turbulence(...)``
         - ``autoPtr<compressibleInterPhaseTransportModel> turbulence_;``
       * - ``fv::options& fvOptions = fv::options::New(mesh);``
         - ``fv::options& fvOptions_;``

   Objects that cannot be default/copy-constructed and are built with
   arguments only known at construction time (e.g. the turbulence model, which
   needs already-built fields) are stored as ``autoPtr<T>`` and allocated in
   the constructor *body* rather than its initializer list (see Step 3).

.. note ::
    You do **not** need to convert every include file's content into a
    class member function. Only the persistent state needs to become member
    data. The imperative code from ``UEqn.H``, ``pEqn.H``, etc. keeps living
    in ``.H`` snippets that get ``#include``-d — just now inside a member
    function body instead of inside ``main()`` (see Step 4).


Step 3 — write the constructor
--------------------------------

The constructor signature is fixed by the run-time selection table:
``compressibleInterFoam(dynamicFvMesh& mesh)``. Its job is to do exactly what
``createTime.H`` + ``createMesh.H`` + ``createFields.H`` used to do, but as a
C++ member-initializer list, in the same order the members were declared in
the header (C++ initializes members in declaration order regardless of the
order written in the list — get this wrong and you silently read
uninitialized data).

.. code :: cpp

    Foam::solvers::compressibleInterFoam::compressibleInterFoam
    (
        dynamicFvMesh& mesh
    )
    :
        solver(mesh),                 // <- always first: constructs mesh_/runTime_
        pimple_(mesh_),
        p_rgh_
        (
            IOobject
            (
                "p_rgh",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        ),
        U_(...),
        phi_(...),
        mixture_(U_, phi_),
        alpha1_(mixture_.alpha1()),
        alpha2_(mixture_.alpha2()),
        rho1_(mixture_.thermo1().rho()),
        rho2_(mixture_.thermo2().rho()),
        rho_(...),
        pMin_("pMin", dimPressure, mixture_),
        g_(...),
        hRef_(...),
        gh_("gh", (g_ & mesh_.C()) + mag(g_)*hRef_),
        ghf_("ghf_", (g_ & mesh_.Cf()) + mag(g_)*hRef_),
        rhoPhi_(...),
        dgdt_("dgdt", alpha1_*fvc::div(phi_)),
        K_("K", 0.5*magSqr(U_)),
        fvOptions_(fv::options::New(mesh_)),
        p_(mixture_.p()),
        T_(mixture_.T()),
        psi1_(mixture_.thermo1().psi()),
        psi2_(mixture_.thermo2().psi()),
        MRF_(mesh_)
    {
        // Anything that needed a "New()" factory call or that depends on
        // several already-built members goes in the constructor body:
        turbulence_.reset
        (
            new compressibleInterPhaseTransportModel
            (
                rho_, U_, phi_, rhoPhi_, alphaPhi10, mixture_
            )
        );
    }

Practical guidance drawn directly from this example:

- ``solver(mesh)`` must be the first entry — it constructs ``mesh_`` and
  ``runTime_``, which every later member depends on.
- Every ``mesh``/``runTime`` reference from the native ``createFields.H``
  becomes ``mesh_``/``runTime_`` (or ``mesh_.time()`` where the native code
  used a bare ``runTime``, since inside the constructor ``runTime()`` the
  accessor is not yet safely callable — use ``mesh_.time()``).
- Reference members (``volScalarField& alpha1_``, ``const volScalarField&
  rho1_``) are bound once here and never reseated afterwards — exactly as
  the native code aliases ``mixture.alpha1()`` into a local reference.
- Anything the native solver built via a factory/``New()`` call that can't sit
  in an initializer list (turbulence model, ``MRF`` in some solvers, etc.)
  is stored as ``autoPtr<T>`` and ``.reset(new T(...))`` in the constructor
  body.


Step 4 — write ``correctPhysics()`` and ``maxDeltaT()``
---------------------------------------------------------

This is the crucial split. The native ``main()`` had two parts:

.. code :: cpp

    while (runTime.run())
    {
        #include "readTimeControls.H"     // <-- becomes maxDeltaT()
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"
        ++runTime;

        while (pimple.loop())             // <-- becomes correctPhysics()
        {
            #include "alphaControls.H"
            #include "compressibleAlphaEqnSubCycle.H"
            turbulence.correctPhasePhi();
            #include "UEqn.H"
            volScalarField divUp(...);
            #include "TEqn.H"
            while (pimple.correct()) { #include "pEqn.H" }
            if (pimple.turbCorr()) turbulence.correct();
        }

        runTime.write();
    }

**``correctPhysics()``** becomes the inner PIMPLE loop, verbatim, with member
names substituted for globals:

.. code :: cpp

    void Foam::solvers::compressibleInterFoam::correctPhysics()
    {
        while (pimple_.loop())
        {
            #include "alphaControls.H"
            #include "compressibleAlphaEqnSubCycle.H"

            turbulence_->correctPhasePhi();

            #include "UEqn.H"
            volScalarField divUp("divUp", fvc::div(fvc::absolute(phi_, U_), p_));
            #include "TEqn.H"

            while (pimple_.correct())
            {
                #include "pEqn.H"
            }

            if (pimple_.turbCorr())
            {
                turbulence_->correct();
            }
        }
    }

**``maxDeltaT()``** overrides the pure virtual from ``solver`` and takes over
the job of ``readTimeControls.H`` / ``CourantNo.H`` / ``alphaCourantNo.H`` /
``setDeltaT.H``: instead of mutating a global ``runTime.setDeltaT(...)``
directly, it *returns* the maximum stable time step, which the outer
GeN-Foam-level ``adjustDeltaT(runTime, solvers)`` collects from every region
solver and applies once, globally:

.. code :: cpp

    Foam::scalar Foam::solvers::compressibleInterFoam::maxDeltaT()
    {
        scalar newDeltaT = mesh_.time().controlDict()
            .lookupOrDefault<scalar>("maxDeltaT", GREAT);

        scalar maxCo = mesh_.time().controlDict().getOrDefault<scalar>("maxCo", 1);

        if (mesh_.time().value() > mesh_.time().controlDict().get<scalar>("deltaT"))
        {
            // ... recompute CoNum / alphaCoNum from phi_, mesh_.V(), etc.
            // ... same arithmetic as native CourantNo.H / alphaCourantNo.H
        }

        return newDeltaT;
    }

Anything from ``main()`` that neither belongs to the per-timestep physics nor
to the time-step estimate (``runTime.write()``, ``runTime.printExecutionTime``,
the outer ``while (runTime.run())`` itself) is **not** ported at all — it is
already handled once, globally, by ``applications/solvers/GeN-Foam/GeN-Foam.C``.
Do not try to reproduce it inside the class.

The remaining virtuals declared in the header
(``correctTightlyCoupledPhysics``, ``correctEnergy``, ``correctFluidMechanics``,
``correctCourant``, ``correctContErr``, ``printContErr``,
``calcCumulContErr``) may simply be left as empty bodies (``{}``) if the
solver has no need for the finer-grained hooks used by some coupling
``loop`` types — this is exactly what ``compressibleInterFoam.C`` does. Only
implement them if a coupling scheme you plan to use actually calls them (grep
``src/multiRegion/loop`` for real callers before spending time on them).

Finally, register the class so ``solver::New()`` can find it by name (this
must be at file scope in the ``.C`` file, right after the includes):

.. code :: cpp

    #include "compressibleInterFoam.H"
    #include "addToRunTimeSelectionTable.H"

    namespace Foam
    {
    namespace solvers
    {
        defineTypeNameAndDebug(compressibleInterFoam, 0);
        addToRunTimeSelectionTable(solver, compressibleInterFoam, dynamicFvMesh);
    }
    }


Step 5 — port the equation files
-----------------------------------

Move every ``.H`` snippet still ``#include``-d from ``correctPhysics()``
(``UEqn.H``, ``pEqn.H``, ``TEqn.H``, ``alphaControls.H``, ``alphaEqn.H``,
``compressibleAlphaEqnSubCycle.H``, ``rhofs.H``, ``alphaSuSp.H``, ...) into
``<solverName>/include/equations/``, and inside each one, rename every
reference to a variable that is now a class member. The rewrite is
mechanical — same physics, same equation, only the identifiers change:

.. list-table:: Typical identifier rewrites inside equation files
    :widths: 40 40
    :header-rows: 1

    * - Native (free variable)
      - Wrapped (member access)
    * - ``rho``
      - ``rho_``
    * - ``U``, ``phi``
      - ``U_``, ``phi_``
    * - ``mesh``
      - ``mesh_`` (or ``mesh()`` accessor)
    * - ``pimple``
      - ``pimple_``
    * - ``turbulence.correct()``
      - ``turbulence_->correct()`` (member is ``autoPtr``, use ``->``)
    * - ``fvOptions``
      - ``fvOptions_``
    * - ``MRF``
      - ``MRF_``

For example, ``UEqn.H`` becomes:

.. code :: cpp

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho_, U_) + fvm::div(rhoPhi_, U_)
      - fvm::Sp(contErr, U_)
      + MRF_.DDt(rho_, U_)
      + turbulence_->divDevRhoReff(U_)
     ==
        fvOptions_(rho_, U_)
    );

    UEqn.relax();
    fvOptions_.constrain(UEqn);

    if (pimple_.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture_.surfaceTensionForce()
                  - ghf_*fvc::snGrad(rho_)
                  - fvc::snGrad(p_rgh_)
                ) * mesh_.magSf()
            )
        );
        fvOptions_.correct(U_);
        K_ = 0.5*magSqr(U_);
    }

Note that local (non-member) variables the equation file declares for its own
scope (``UEqn`` itself, ``rAU``, ``phiHbyA`` in ``pEqn.H``, etc.) are left
completely untouched — they stay ordinary stack locals inside
``correctPhysics()``, exactly like in the native solver, since they don't need
to persist across calls.


Step 6 — wire up the build (``Make/files`` and ``Make/options``)
--------------------------------------------------------------------

The single biggest build-system difference: the native solver builds an
**executable** (``EXE``); the wrapped class builds into a **shared library**
(``LIB``) that GeN-Foam/OFFBEAT dynamically link against and select from at
run time.

Native ``compressibleInterFoamNATIVE/Make/files``:

.. code :: text

    compressibleInterFoam.C

    EXE = $(FOAM_APPBIN)/compressibleInterFoam

Wrapped ``applications/modules/openFoamImported/Make/files`` (shared across
every solver class placed under ``openFoamImported/``):

.. code :: text

    compressibleInterFoam/compressibleInterFoam.C

    LIB = $(FOAM_USER_LIBBIN)/libOFSolvers

Add one line per new solver's ``.C`` file to this shared ``Make/files`` as you
add solvers to the module directory.

``Make/options`` needs, in addition to whatever native include paths /
libraries the original solver required:

- ``-I./include`` and ``-I./<solverName>/include/equations`` so the class
  source can find its own headers and equation snippets,
- ``-I<repo>/src/multiRegion/lnInclude`` so it can see ``solver.H``,
- any auxiliary OpenFOAM libraries the native solver needed, referenced either
  through the standard ``$(LIB_SRC)`` (for stock OpenFOAM libraries) or through
  ``$(FOAM_SOLVERS)/<category>/<solverName>/<lib>/lnInclude`` when the native
  solver ships its own private library (as ``compressibleInterFoam`` does for
  ``twoPhaseMixtureThermo`` and ``VoFphaseCompressibleTurbulenceModels`` —
  these are *not* recompiled inside this repository; the module reuses the
  copies already built as part of the OpenFOAM installation).

.. code :: text

    EXE_INC = \
        -I./.. \
        -I./include \
        -I./compressibleInterFoam \
        -I./compressibleInterFoam/include/equations \
        -I./../../../src/multiRegion/lnInclude \
        -I./../../../src/fvPatchFields/lnInclude \
        -I$(FOAM_SOLVERS)/multiphase/compressibleInterFoam/twoPhaseMixtureThermo/lnInclude \
        -I$(FOAM_SOLVERS)/multiphase/compressibleInterFoam/VoFphaseCompressibleTurbulenceModels/lnInclude \
        ... (standard finiteVolume/thermophysicalModels/turbulenceModels includes)

    LIB_LIBS = \
        -L$(FOAM_USER_LIBBIN) \
        -lfiniteVolume -lfvOptions -lmeshTools \
        -ltwoPhaseMixtureThermo -ltwoPhaseSurfaceTension \
        ...

.. note ::
    The ``compressibleInterFoamNATIVE`` directory is deliberately left with
    its own, separate ``Make/files``/``Make/options`` pointing at
    ``EXE = $(FOAM_APPBIN)/...`` — but it is never invoked by any top-level
    ``Allwmake``. It exists purely so a developer can still ``wmake`` the
    stock solver standalone to sanity-check behaviour against the wrapped
    version if ever in doubt.


Step 7 — select the solver from a case
-----------------------------------------

Once ``libOFSolvers`` is rebuilt, the class is available to
``solver::New()`` by its ``TypeName``. Select it for a region in
*system/regionsDict* exactly like any other foamForNuclear solver
(see :ref:`Achieving coupled solutions <couplingGF>`):

.. code :: cpp

    // In system/regionsDict

    regionSolvers
    {
        Level_0
        {
            fluidRegion    compressibleInterFoam;
        }
    }

Optionally, add a ``<solverName>.yaml`` file next to the class (see
``compressibleInterFoam.yaml``) describing it for the auto-generated API
documentation:

.. code :: yaml

    description: |-
      Short description of the solver, shown in the generated docs.

    admonitions:
      - kind: warning
        body: |-
          Preliminary autogenerated documentation. Please refer to the
          source code and tutorials for more reliable information.

    options:
      - key: <parameter>
        type: word
        required: true
        default:
        description: <description>


Step 8 — build and test
--------------------------

Rebuild the module (``wmake`` from ``applications/modules/openFoamImported``,
or the repository's top-level ``Allwmake``), then run any tutorial case whose
``regionsDict`` selects your new solver name and confirm it starts, iterates,
and writes fields exactly as the native application did on the same case
(comparing residuals/field values against a run of the untouched
``compressibleInterFoamNATIVE`` copy, built standalone with plain ``wmake``,
is the most reliable regression check while porting).


Summary — concept-to-concept mapping
---------------------------------------

.. list-table::
    :widths: 45 55
    :header-rows: 1

    * - Native OpenFOAM solver concept
      - foamForNuclear class equivalent
    * - ``main(argc, argv)``
      - constructor + ``correctPhysics()`` + ``maxDeltaT()`` (no ``main`` of
        its own — driven by ``GeN-Foam.C``)
    * - ``createTime.H`` / ``createMesh.H``
      - handled by ``solver::solver(dynamicFvMesh&)``, called first in the
        initializer list
    * - ``createFields.H`` (and other ``create*.H``)
      - private member variables (trailing ``_``), built in the constructor's
        initializer list / body
    * - global fields/objects (``U``, ``rho``, ``turbulence``, ``mesh``, ...)
      - member variables (``U_``, ``rho_``, ``turbulence_``, ``mesh_``)
    * - ``while (pimple.loop()) { ... }`` body
      - ``correctPhysics()`` override
    * - ``readTimeControls.H`` / ``CourantNo.H`` / ``setDeltaT.H``
      - ``maxDeltaT()`` override (returns a value instead of mutating
        ``runTime`` directly)
    * - ``EXE = $(FOAM_APPBIN)/solverName``
      - ``LIB = $(FOAM_USER_LIBBIN)/libOFSolvers``
    * - selecting the solver: running the executable
      - selecting the solver: ``regionSolvers/Level_0`` entry in
        *system/regionsDict*, resolved via ``solver::New()`` and
        ``addToRunTimeSelectionTable``
