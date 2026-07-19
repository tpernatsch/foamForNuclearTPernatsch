Post-processing
===============

.. warning::

   This page is a work in progress.

-------------------------------------------------------------------------------

Graphical visualisation using ParaView
--------------------------------------

Several tutorials exist on how to use ParaView for visualising CFD
simulation results. This guide assumes that the reader already has a
basic familiarity with ParaView and with post-processing OpenFOAM
simulations.

A fairly comprehensive overview is provided in the
`OpenFOAM User Guide <https://cfd.direct/openfoam/user-guide/v9-paraview/>`_.

The recommended approach for visualising the solution fields from
OFFBEAT simulations is to use ``paraFoam``, the ParaView plugin that is
shipped with OpenFOAM. Although ParaView includes its own built-in
OpenFOAM reader, ``paraFoam`` generally provides better support for
cases containing ``cellZones`` (OFFBEAT uses ``cellZones`` to define
different material regions in the mesh).

Using ``paraFoam`` makes it possible to focus on specific materials
(e.g. fuel or cladding), or to visualise quantities such as the gap
width or the fuel-cladding contact pressure.

The user launches ParaView from the command line as follows:

.. code-block:: bash

   paraFoam -case <case_directory>

In some cases ``paraFoam`` may not be available, for example when running
simulations on a computing cluster. In that case the ParaView built-in
OpenFOAM reader can be used:

* create an empty file named ``.foam`` in the case folder
* open this file directly in ParaView.

-------------------------------------------------------------------------------

List of fields
--------------

The list of fields available for visualisation using ParaView is quite
extensive. The following provides a brief overview.

Basic solution fields
^^^^^^^^^^^^^^^^^^^^^

* displacement vector ``D`` [m]
* temperature ``T`` [K]

Strain
^^^^^^

* total strain tensor ``epsilon`` [-]
* creep strain ``epsilonCreep`` [-]
* additional strain components depending on the model

Stress
^^^^^^

* stress tensor ``sigma`` [Pa]
* von Mises stress ``sigmaEq`` [Pa]
* hydrostatic stress ``sigmaHyd`` [Pa]
* yield stress ``sigmaY`` [Pa]

Other quantities
^^^^^^^^^^^^^^^^

* burnup ``BU`` [MWd/MTU]
* fast fluence ``fastFluence``

Gap-related quantities
^^^^^^^^^^^^^^^^^^^^^^

* gap width ``gapWidth`` [m]
* gap contact pressure ``interfaceP`` [Pa]

These quantities are typically non-zero only on the boundaries of the
gap.

-------------------------------------------------------------------------------

TUBRNP-related fields
---------------------

Examples include:

* isotope concentrations (e.g. ``N_Am241``, ``N_Pu240``)
* burnup in MWd/kgU (``BUkgU``)
* radial form factor (``formFactorTUBRNP``)

-------------------------------------------------------------------------------

SCIANTIX-related fields
-----------------------

Examples include:

* ``Grain_radius``
* ``Gas_produced``
* ``Gas_grain``
* ``Gas_grain_solution``
* ``Gas_grain_bubbles``
* ``Gas_boundary``
* ``Gas_released``
* ``Effective_burn_up``
* ``Oxygen_to_metal_ratio``
* ``Helium_produced``
* ``Helium_grain``
* ``Helium_grain_solution``
* ``Helium_grain_bubbles``
* ``Helium_boundary``
* ``Helium_released``
* ``Intragranular_bubble_concentration``
* ``Intragranular_bubble_radius``
* ``Intragranular_gas_swelling`` (tensor)
* ``Intergranular_bubble_concentration``
* ``Intergranular_atoms_per_bubble``
* ``Intergranular_vacancies_per_bubble``
* ``Intergranular_bubble_radius``
* ``Intergranular_bubble_area``
* ``Intergranular_bubble_volume``
* ``Intergranular_fractional_coverage``
* ``Intergranular_saturation_fractional_coverage``
* ``Intergranular_gas_swelling`` (tensor)
* ``Intergranular_fractional_intactness``

-------------------------------------------------------------------------------

Tips and tricks for graphical post-processing
---------------------------------------------

Dealing with the high aspect ratio of fuel rod models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fuel rod models typically have a very high aspect ratio
(height ≫ radius). As a consequence, results can appear as extremely
thin vertical structures in ParaView and radial behaviour becomes
difficult to visualise.

In such cases the **Transform** filter can be used to scale the
z-direction.

For typical LWR fuel rods with approximately

.. math::

   \frac{R}{H} \approx 10^{-3}

a scaling of approximately

.. math::

   z \approx 5 \times 10^{-3}

generally provides a good visualisation.

The option **Transform all Input Vectors** should be disabled so that
vector and tensor fields are not distorted.

If the user wishes to visualise displacement and gap closure, the
**Warp By Vector** filter can be applied using the displacement field
``D`` (or ``DD``).

* warp factor = 1 → realistic deformation
* larger factors → enhanced visualisation of deformation

This filter should be applied **before** the Transform filter to avoid
unrealistic distortion in the axial direction.

-------------------------------------------------------------------------------

Saving ParaView states
----------------------

ParaView state files can be used to save the full configuration of a
visualisation pipeline.

This avoids having to reopen the case and recreate filters every time
ParaView is launched.

See the ParaView documentation:

https://www.paraview.org/Wiki/Advanced_State_Management#Save_State_/_Load_State

-------------------------------------------------------------------------------

Post-processing using OpenFOAM utilities
----------------------------------------

Graphical visualisation with ParaView is useful for qualitative
analysis, but quantitative analysis is usually performed using
OpenFOAM sampling utilities.

Descriptions of these utilities are available in the OpenFOAM User
Guide:

* https://cfd.direct/openfoam/user-guide/v9-post-processing-cli
* https://cfd.direct/openfoam/user-guide/v9-graphs-monitoring

Two utilities are particularly relevant for fuel behaviour applications.

-------------------------------------------------------------------------------

Probing data
------------

Probes can be defined in ``controlDict`` to extract time-dependent
values at specific mesh locations.

These are useful in validation exercises where experimental
measurements (e.g. thermocouples) are available.

Example: probing centreline temperature and burnup at three axial
locations.

``system/probes``

.. code-block:: c

   points
   (
       (0.0 0.0 1.0)
       (0.0 0.0 1.5)
       (0.0 0.0 2.0)
   );

   fields (T BU);

   #includeEtc "caseDicts/postProcessing/probes/probes.cfg"

``system/controlDict``

.. code-block:: c

   functions
   {
       #includeFunc probes
   }

During the simulation, OFFBEAT writes ASCII tables in

``postProcessing/probes/<time>/``

These files can easily be plotted using any external tool.

-------------------------------------------------------------------------------

Line graphs
-----------

Radial and axial profiles (e.g. temperature ``T``) can be extracted
using OpenFOAM sampling tools.

Example configuration in ``system/sampleLines``:

.. code-block:: c

   type                sets;
   interpolationScheme cellPointFace;
   setFormat           csv;
   libs                ("libsampling.so");

   sets
   (
       radialProfile
       {
           type    lineCellFace;
           axis    distance;
           start   (0 0 1.6);
           end     (5.6e-3 0 1.6);
       }

       axialProfile
       {
           type    lineCellFace;
           axis    distance;
           start   (0 0 0);
           end     (0 0 3.2);
       }
   );

   fields (T);

After the simulation finishes:

.. code-block:: bash

   postProcess -func sampleLines

This generates CSV files in

``postProcessing/sampleLines/<time>/``

which can be plotted using standard graphing tools.