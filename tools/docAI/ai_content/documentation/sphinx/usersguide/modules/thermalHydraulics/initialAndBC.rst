.. _userguide_thermalhydraulics_initialAndBC:


Initial and boundary conditions
-------------------------------

Initial and boundary conditions follow the usual OpenFOAM logic for one- and two-phase solvers. OpenFOAM provides most of the boundary conditions needed for thermal-hydraulics models. In addition, foamForNuclear includes a few boundary conditions in its thermal-hydraulics module:

.. toctree::
   :maxdepth: 1
   :glob:

   ../../../cppapi/generated/fvPatchFields/thermalHydraulics/**


.. raw:: html

   <br><br>


Power densities can be provided directly to the :ref:`onePhase <onePhase>`, :ref:`onePhaseLegacy <onePhaseLegacy>`, and :ref:`twoPhase <twoPhase>` solvers as fields in the initial time folder. There are two types of power densities:

   - ``powerDensityLiquid`` is given to the liquid.
   - ``powerDensityStructure`` is given to the structure and used by the :doc:`../../../../../cppapi/generated/porousMediaModels/phaseModels/structureModels/powerModels/powerModel` (see :ref:`structureProperties  <modules_thermalHydraulics_porousMedium_structureProperties>`).

.. note::

    In two-phase simulations with liquid fuel, the ``powerDensity`` is applied to anything that is liquid in the thermal-hydraulics model. You are supposed to have one liquid phase and one gas phase. Otherwise, the power will be counted twice.


.. note::

    The power densities in the thermal-hydraulics sub-solver are ALWAYS the physical ones. For instance, when the *nuclearFuelPin* model is used for pin-based reactors, *powerDensity* refers to the power density inside the fuel matrix. For liquid fuel, *powerDensity* is the power density in the liquid. They are not the power densities smeared over the whole volume.


.. note::

    If you calculate the ``powerDensity`` using, for example, Serpent, you have to divide it by the fuel fraction before feeding it to GeN-Foam (see :ref:`Neutronics <neutronics_the-reactorstate-dictionary>`).