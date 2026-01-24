.. _userguide_thermalhydraulics_settingcase:

Setting a case
--------------

Initial and boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Initial and boundary conditions adopt the usual OpenFOAM logic for one- and two-phase solvers. 
OpenFOAM provides most of the boundary conditions one may need for
thermal-hydraulics models. In addition, a few boundary conditions have been
included in the thermal-hydraulics module of foamForNuclear:

toctreeHere

powerDensities can be provided to the :ref:`onePhase <onePhase>`, :ref:`onePhaseLegacy <onePhaseLegacy>` and :ref:`twoPhase <twoPhase>` solvers direclty as fields in the initial time folder. There two types of power densities:

   - ``powerDensity`` is given to the liquid
   - ``powerDensityStructure`` is given to the structure and used by the :doc:`../../../../../cppapi/generated/porousMediaModels/phaseModels/structureModels/powerModels/powerModel` (see :ref: `structureProperties  <modules_thermalHydraulics_porousMedium_structureProperties>`. 

.. note::

    In two-phase simulations with liquid fuel, the powerDensity to anything that is liquid in thermal-hydraulics. You are supposed to
    have one liquid and one gas. Otherwise, power will be counted twice.


.. note::

    The power densities in the thermal-hydraulic sub-solver are ALWAYS the
    physical ones: for instance, when the *nuclearFuelPin* model is used for
    pin-based reactors, *powerDensity* refers to the power density inside the
    fuel matrix. For liquid fuel, the *powerDensity* is the power density in the
    liquid. They are not the power densities smeared over the whole volume.

.. note::

    If  calculating the powerDensity using e.g. Serpent, one has to divide it by
    the fuel fraction before feeding it to GeN-Foam (see
    :ref:`Neutronics <neutronics_the-reactorstate-dictionary>`)


Discretization and solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Details for discretization and solution of equations are handled in a standard
OpenFOAM way, i.e., through the *fvSolution* and *fvSchemes* dictionaries in
*system/[nameOfTheFluidRegion]*.

