.. _userguide_neutronics_nuclearData:


The *nuclearData* dictionary
----------------------------

In GeN-Foam, cross-sections and several other neutronics properties are handled
by the :ref:`XS.H <XS>` class. Detailed explanations on the file format are provided in
:ref:`XS.H <XS>` and in the tutorials (e.g `3D_SmallESFR
<https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/neutronicsProperties>`_).

The *nuclearData* dictionary can be found under *constant/(neutronicsRegionName)/neutronicsProperties*. It
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
<https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/neutronicsProperties>`_).

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

.. figure:: ../../../images/fig_RBF_interpolationOrders.png
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

    nuclearData
    {
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
    }


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

    nuclearData
    {
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
    }

One can find more details on all the parameters in the :ref:`XS.H <XS>` file and commented
examples of *nuclearData* in the tutorials
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/neutronicsProperties>`_ (for diffusion or SP3),
`Godiva_SN <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/Godiva_SN/constant/neutroRegion/neutronicsProperties>`_ (for discrete ordinates) and
`2D_onePhaseAndPointKineticsCoupling <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/2D_onePhaseAndPointKineticsCoupling/rootCase/constant/neutroRegion/neutronicsProperties>`_ (for point kinetics).
`2D_onePhaseAndSubcriticalPointKineticsCoupling <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/featureCases/2D_onePhaseAndSubcriticalPointKineticsCoupling/rootCase/constant/neutroRegion/externalSource>`_ (for subcritical point kinetics).