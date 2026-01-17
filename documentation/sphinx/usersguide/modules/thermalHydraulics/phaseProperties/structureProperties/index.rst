.. _modules_thermalHydraulics_porousMedium_regimeMapModels:

The *structureProperties* sub-dictionary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A coarse-mesh porous-medium treatment  implies that the cell zone is
modeled without resolving the sub-scale structure (e.g., the fuel rods or the
heat exchanger tubes). This means that dedicated models are needed to describe their behavior. foamForNuclear allows
modeling simultaneously in the same region both a "passive
structure" and a "power model" . Passive
structures are simply modeled as a heat capacity and they can be used for
instance to model assembly wrappers or reflector structures. In essence, the
fluid will interact cell-by-cell with these passive structures: they will take
energy from the fluid if the temperature of the fluid is higher than that of the
surface of the structure, and vice versa. Power models are more compelex and are typically used to model components that produce power (either nuclear or electrical).  An example of a power model is the
nuclearFuelPin, which can be used to model a standard pin-type fuel. This power
model is capable of getting the power density from neutronics, solving a 1-D
model for heat transfer in the fuel, and giving back to the fluid the
temperature at the surface of the cladding. The fluid will then be capable of
calculating the heat transfer with the fuel based on the cladding surface
temperature and the Nusselt number. 

The  *structureProperties* dictionary consists of a series of sub-dictionaries
whose properties are applied to the cellZones that have their name
in the subDict keys. The presence of any of these
sub-dictionaries is not mandatory, and the structure defaults, globally, to
having a null volumeFraction and infinite (1e6m) hydraulic diameter.
For convenience, if multiple cellZones share the
same exact properties, these names can be grouped by using a colon as
separator. Two examples of *structureProperties* dictionaries are reported below. *diagrid*, *axialReflector*, etc. are name of cell zones. The main run-time selectable models are:

.. toctree::
   :maxdepth: 1

   powerModels
   powerOffCriterionModels


.. code :: cpp

    structureProperties
    {
        "diagrid:axialReflector:radialReflector:follower:controlRod"
        {
            volumeFraction      0.718520968;
            Dh                  0.00365; //hydraulic diameter

            passiveProperties 
            {
                volumetricArea  5;
                rhoCp           4.8e6;
                T               668;
            }

        }

        "innerCore"
        {
            volumeFraction      0.718520968;
            Dh                  0.00365;

            powerModel 
            {
                type                nuclearFuelPin;

                powerDensity        0;  //- fields on disk have priority, if they
                                        //  are not found, this value is used
                fuelInnerRadius     0.0012;
                fuelOuterRadius     0.004715;
                cladInnerRadius     0.004865;
                cladOuterRadius     0.005365;
                fuelMeshSize        30;
                cladMeshSize        5;
                fuelRho             10480;
                fuelCp              250;
                cladRho             7500;
                cladCp              500;
                gapH                3000;
                fuelK               3;
                cladK               20;
                fuelT               668;
                cladT               668;
            }

            passiveProperties // these are the properties of the metallic wrappers
            {
                volumetricArea  5;
                rhoCp           4.8e6;
                T               668;
            }

        }
    }


.. code :: cpp

    structureProperties
    {
        "lowIn:topIn"
        {
            volumeFraction  0.523128;
            Dh              0.005469;
            localTortuosity (0.6366197724 0.6366197724 1.0);
            localDhAnisotropy     (1 1 1);
            localX       (1 0 0);
            localZ       (0 0 1);

            powerModel
            {
                type            heatedPin;
                innerRadius     0;
                outerRadius     0.003;
                k               20;
                meshSize        6;
                T               653.15;
                Cp              500;
                rho             7700;

                powerOffCriterionModel
                {
                    type            timeThreshold;
                    time            12.475;
                }
            }
        }

    }






