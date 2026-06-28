Gap gas
=======

The ``gapGas`` model handles the evolution of gap volume, pressure,
temperature, and gas composition in the thermo-mechanics module
(OFFBEAT / GeN-Foam solid region).

-------------------------------------------------------------------------------

.. rubric:: Selection in ``solverDict``

The gap gas model is selected in ``solverDict`` via the ``gapGas``
dictionary (located in ``constant/`` for OFFBEAT, or in
``constant/<regionName>/`` for GeN-Foam):

.. code-block:: cpp

   gapGas
   {
       type Frapcon;
       // additional options (if any)
   }

If the ``gapGas`` dictionary is not present in ``solverDict``, the gap
gas behaviour is neglected.

-------------------------------------------------------------------------------

.. warning::

   In the ``gapGas`` model, the gap volumes refer to the fraction of the
   rod or particle actually represented in the mesh.

   For example, consider a rod plenum volume of ``36 cm^3``. If the rod
   is represented using a ``2°`` wedge mesh, the plenum volume represented
   in the simulation becomes

   .. math::

      V_{\text{model plenum}} =
      \frac{2}{360} \cdot \frac{36}{10^6} \, m^3
      = 2 \times 10^{-7} \, m^3

   For TRISO fuels, consider a spherical gap of ``36 cm^3``. If the sphere
   is represented using a ``2°`` wedge mesh (1D representation), the fraction
   of the sphere represented in the model is

   .. math::

      \text{fraction} =
      \frac{(2\pi/180)^2}{4\pi}
      = 9.696273622 \times 10^{-5}

   The corresponding gap volume in the TRISO model is therefore

   .. math::

      V_{\text{model}} =
      36 \times 10^{-6}
      \cdot
      9.696273622 \times 10^{-5}
      \, m^3
      =
      3.49 \times 10^{-9} \, m^3

-------------------------------------------------------------------------------

The thermo-mechanics module currently provides the following gap gas models:

.. toctree::
   :maxdepth: 1
   :caption: Gap gas models
   :glob:

   ../../../../../cppapi/generated/offbeatLib/gapGasModel/*