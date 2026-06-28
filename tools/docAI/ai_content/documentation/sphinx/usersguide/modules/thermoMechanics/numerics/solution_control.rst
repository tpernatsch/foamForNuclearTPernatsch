Solution control (``fvSolution``)
=================================

The ``fvSolution`` file defines how OFFBEAT handles the linear and nonlinear
solution of its governing equations. Although OFFBEAT builds on the standard
OpenFOAM infrastructure, this document provides an integrated description that
does not assume full familiarity with OpenFOAM.

The goal is to provide enough background to understand the structure of the
file, the purpose of each entry, and the practical considerations behind
typical parameter choices.

.. note::

   Users who need exhaustive details on solver algorithms, preconditioners,
   or advanced linear-algebra options may refer to the official OpenFOAM
   documentation.

-------------------------------------------------------------------------------

Overview of the file structure
---------------------------------

A typical ``fvSolution`` file for OFFBEAT contains three main sections.

``solvers { … }``
^^^^^^^^^^^^^^^^^

This block defines how the linear system arising from the discretised
equation should be solved for each field.

Each entry specifies:

* the solver type (e.g. ``PCG``, ``PBiCGStab``, ``GAMG``)
* the preconditioner
* the linear convergence criteria

In a multiphysics simulation, each field (temperature, displacement, fluxes,
etc.) produces its own linear system and therefore has its own solver entry.

Nonlinear solution controls (``stressAnalysis { … }``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While the ``solvers`` block governs how *linear* systems are solved, the
physics block defines how the *nonlinear* system is handled.

OFFBEAT uses a segregated Picard-type strategy:

* fields are solved sequentially
* the procedure repeats until convergence criteria are satisfied

This block controls:

* number of correctors
* maximum number of outer loops
* nonlinear convergence criteria

``relaxationFactors { … }``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Under-relaxation is a standard technique for stabilising segregated
nonlinear algorithms.

Relaxation may be applied to:

* **fields** (solution values)
* **equations** (system matrix)

Mechanically coupled problems may require additional damping, whereas
purely transient thermal simulations may require little or no relaxation.

Together, these blocks determine both the numerical robustness and the
computational performance of a simulation.

-------------------------------------------------------------------------------

``solvers {}`` — linear algebra settings
-------------------------------------------

Every discretised equation in OFFBEAT results in a linear system

.. math::

   A x = b

The ``solvers`` dictionary specifies how this system is solved.

Each field has its own entry (e.g. ``D``, ``T``).

What each entry contains
^^^^^^^^^^^^^^^^^^^^^^^^

``solver``

The iterative method used to solve the linear system.

* ``PCG`` – appropriate for symmetric positive-definite systems
* ``PBiCGStab`` – suitable for general non-symmetric systems
* ``GAMG`` – algebraic multigrid solver for large or stiff problems

``preconditioner``

An approximation to :math:`A^{-1}` used to accelerate convergence.

* ``DIC`` / ``FDIC`` – incomplete Cholesky for symmetric systems
* ``DILU`` – incomplete LU for general systems
* ``GAMG`` uses internal smoothers instead

``tolerance``

Absolute residual threshold for stopping the linear solver.

``relTol``

Relative stopping criterion based on reduction of the residual.

Setting ``relTol 0;`` disables this criterion.

``minIter`` / ``maxIter``

Minimum and maximum number of solver iterations.

.. note::

   The residual reported by the linear solver is **not** the same as the
   nonlinear residual used in the ``stressAnalysis`` block.

   * Linear residual → accuracy of the linear solve
   * Nonlinear residual → consistency between Picard iterations

Typical parameter ranges
^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1

   * - System type
     - Solver
     - Preconditioner
     - ``tolerance``
     - ``relTol``
     - Comments
   * - Symmetric diffusion / elasticity
     - ``PCG``
     - ``DIC`` / ``FDIC``
     - ``1e-8 – 1e-10``
     - ``0`` or ``1e-3``
     - Common in thermo-mechanics
   * - Non-symmetric systems
     - ``PBiCGStab``
     - ``DILU``
     - ``1e-8 – 1e-10``
     - ``0 – 1e-3``
     - More robust for mixed systems
   * - Very large systems
     - ``GAMG``
     - internal
     - ``1e-8 – 1e-10``
     - ``0 – 1e-3``
     - Best performance when applicable

-------------------------------------------------------------------------------

``stressAnalysis {}`` — nonlinear solution controls
------------------------------------------------------

The nonlinear coupling between fields is handled by a Picard
(segregated) loop.

At each outer iteration:

1. each field is solved sequentially
2. the most recent values of other fields are used
3. residuals are evaluated

-------------------------------------------------------------------------------

How nonlinear residuals work
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After each corrector, OFFBEAT evaluates nonlinear residuals.

* **absolute residual**

  Standard OpenFOAM residual measuring equation imbalance.

* **relative residual**

  Change between successive iterations relative to the previous
  time-step.

In addition, optional **absolute error criteria** may be used:

* ``absErrD``
* ``absErrT``

These measure the difference between two successive iterations.

Convergence logic
"""""""""""""""""

By default, the outer loop terminates when both:

* absolute residual criteria
* relative residual criteria

are satisfied for all monitored fields.

-------------------------------------------------------------------------------

Reference constraints
^^^^^^^^^^^^^^^^^^^^^

In some cases, it is useful to fix the value of a variable at selected
cells.

OFFBEAT provides:

* ``referencePairs`` – fixes **displacement**
* ``referencePairsThermal`` – fixes **temperature**

Each list contains tuples:

* ``(cellID, vector)`` for displacement
* ``(cellID, scalar)`` for temperature

Example

.. code-block:: c

   stressAnalysis
   {
       referencePairs          ( (12345 (0 0 0)) (67890 (0 0 0)) );
       referencePairsThermal   ( (13579 600) );
   }

-------------------------------------------------------------------------------

``nCorrectors``
^^^^^^^^^^^^^^^

Controls the number of Picard correctors per outer iteration.

Scalar form:

.. code-block:: c

   nCorrectors 3;

Dictionary form:

.. code-block:: c

   nCorrectors
   {
       default 1;
       D       6;
       T       2;
   }

Typical values range from **3–10**. Larger values improve coupling at the
cost of additional linear solves.

-------------------------------------------------------------------------------

Summary of keys in ``stressAnalysis``
-------------------------------------

.. list-table::
   :header-rows: 1

   * - Key
     - Meaning
   * - ``nCorrectors``
     - Number of Picard correctors per outer iteration
   * - ``maxOuterIter``
     - Maximum number of outer iterations
   * - ``<field>``
     - Absolute nonlinear residual target
   * - ``rel<field>``
     - Relative nonlinear residual
   * - ``useRelRes<field>``
     - Enables or disables relative residual criterion
   * - ``absErr<field>``
     - Absolute error threshold
   * - ``referencePairs``
     - Fixed displacement locations
   * - ``referencePairsThermal``
     - Fixed temperature locations

-------------------------------------------------------------------------------

``relaxationFactors {}`` — stabilisation mechanisms
------------------------------------------------------

Segregated algorithms may exhibit oscillations during nonlinear iterations.

Under-relaxation helps stabilise the solution.

Field relaxation
^^^^^^^^^^^^^^^^^

After solving the linear system, the solution is blended:

.. math::

   x^{n+1} = \alpha x_\text{new} + (1-\alpha)x^n

Lower values of :math:`\alpha` increase stability but slow convergence.

Equation relaxation
^^^^^^^^^^^^^^^^^^^

This relaxes the assembled matrix coefficients before the linear solve.

It is useful when the system is nearly singular, for example in:

* poorly constrained mechanical problems
* pellet stacks
* 2D disc simulations

Typical values: **0.9–0.99**.

-------------------------------------------------------------------------------

Annotated example
------------------

.. code-block:: c

   solvers
   {
       "D|DD"
       {
           solver          PCG;
           preconditioner  FDIC;
           tolerance       1e-10;
           relTol          1e-3;
       }

       T
       {
           solver          PCG;
           preconditioner  FDIC;
           tolerance       1e-10;
           relTol          1e-3;
       }
   }

   stressAnalysis
   {
       nCorrectors
       {
           default 1;
           D       6;
       }

       maxOuterIter 1000;

       referencePairs        ();
       referencePairsThermal ();

       D  (1e-5 1 1e-5);
       T  1e-5;
   }

   relaxationFactors
   {
       fields
       {
           D 0.9;
       }

       equations
       {
           // D 0.99;
       }
   }