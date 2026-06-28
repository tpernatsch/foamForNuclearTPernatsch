.. _userguide_neutronics_crmove:

The *CRMove* dictionary
-----------------------

The *CRMove* dictionary is located under *constant/neutroRegion/*. It is optional and provides input data for control rod movement in spatial neutronics sub-solvers. Control rods can be moved from their initial position to a new one by specifying the initial and final times of insertion/extraction, as well as the insertion/extraction speed (positive speed for insertion).

A commented example is available in the tutorial
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/CRmove>`_. However, this option is not actually used in the tutorial.

The following keywords can be used to define control rod movement:

:startTime: Time at which the control rod should start moving.
:endTime: Time at which the control rod should stop moving.
:speed: Speed of the control rod movement in m/s.
:initialDistanceFromMeshCR: Initial distance of the control rod relative to the initial position defined in the mesh as a cellZone. A positive value means the control rod is offset along the -Z axis.
:followerName: Name of the cellZone representing the domain that replaces the control rod during the movement out of its initial cellZone in the mesh.

.. code :: cpp

    zones
    (
        controlRod1
        {
            startTime       1001;
            endTime         1003;
            speed           0.5; // positive entering the core
            initialDistanceFromMeshCR 0.0; // distance from the position defined in the mesh (postive towards the core)
            followerName    follower; // name of the cellZone representing the follower
        }

        controlRod2
        ...
    );