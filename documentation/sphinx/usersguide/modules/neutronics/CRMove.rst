.. _userguide_neutronics_crmove:

The *CRMove* dictionary
-----------------------

The *CRMove* dictionary can be found under *constant/neutroRegion/*. It is optinal 
and contains input data for control rod movement in spatial neutronics solvers. 
Control rods can be moved from the initial position to a new one by selecting the 
initial and final time of the insertion/extraction and the speed of 
insertion/extraction (positive speed for insertion).

One can find a commented example in the tutorial
`3D_SmallESFR <https://gitlab.com/foamForNuclear/foamForNuclear/-/tree/master/tutorials/reactorCases/3D_SmallESFR/extendedThermoMechanics/constant/neutroRegion/CRmove>`_, though this option is not actually used in the tutorial.

Here is the list of keywords that can be used to set a control rod movement:

:startTime: Time when the control rod should start moving.
:endTime: Time when the control rod should stop moving.
:speed: Speed of the control rod movement in m/s.
:initialDistanceFromMeshCR: Initial distance of control rod with respect to the
                            initial position defined in the mesh as a cellZone.
                            A positive value means the control rod is offset
                            along the -Z axis.
:followerName: Name of the cellZone representing the domain that replaces the
               control rod during the movement out of its initial cellZone in
               the mesh.

.. code :: cpp

    zones
    (
        controlRod1
        {
            startTime       1001;
            endTime         1003;
            speed           0; // positive entering the core
            initialDistanceFromMeshCR 0.0; // distance from the position defined in the mesh (postive towards the core)
            followerName    follower; // name of the cellZone representing the follower
        }

        controlRod2
        ...
    );
