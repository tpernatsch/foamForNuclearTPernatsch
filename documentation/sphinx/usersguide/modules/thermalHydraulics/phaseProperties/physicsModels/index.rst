.. _modules_thermalHydraulics_porousMedium_physicsModels:

The *physicsModels* sub-dictionary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The *physicsModels* sub-dictionary is used in one-phase simulations to define the models that describe the effect of the structure on the fluid flow.  Two types of models can be defined: drag models via the ``dragModels`` sub-dictionary and heat transfer models via the ``heatTransferModels`` sub-dictionary. Each of these two sub-dictionaries consists of a series of other sub-dictionaries whose properties are applied to the cell zones that have their name in the sub-dictionary keys. For both the drag models and the heat transfer models, three options exist:

- Use ``constant`` as  type and give a constant value using the keyword ``value``.
- Use as type a single model chosen among those available (:ref:`drag models <modules_thermalHydraulics_porousMdeium_FSdrag>`, :ref:`heat transfer models <modules_thermalHydraulics_porousMedium_FSHeatTransferCoefficientModels>`)
- Use ``byRegime`` as  type, and then chose one of the available models (:ref:`drag models <modules_thermalHydraulics_porousMdeium_FSdrag>`, :ref:`heat transfer models <modules_thermalHydraulics_porousMedium_FSHeatTransferCoefficientModels>`) for every regime. This requires having set a :ref:`regime map <modules_thermalHydraulics_porousMedium_regimeMapModels>`.



Two examples of dictionaries are reported below. 

.. code :: cpp
   
   physicsModels
   {
      dragModels
      {
         "diagrid:axialReflector:radialReflector:follower:controlRod:innerCore:outerCore"
         {
               type    ReynoldsPower;
               coeff   0.687;
               exp     -0.25;
         }
      }

      heatTransferModels
      {
         "diagrid:axialReflector:radialReflector:follower:controlRod:innerCore:outerCore"
         {
               type        byRegime;
               regimeMap   "lamTurb";

               "laminar"
               {
                  // Nu = const + coeff * Re^expRe * Pr^expPr
                  type    NusseltReynoldsPrandtlPower;
                  const   4;
                  coeff   0;
                  expRe   0;
                  expPr   0;
               }
               "turbulent"
               {
                  type    NusseltReynoldsPrandtlPower;
                  const   4.82;
                  coeff   0.0185;
                  expRe   0.827;
                  expPr   0.827;
               }
         }
      }
   }


.. code :: cpp

    physicsModels
    {
        dragModels
        {
            [nameOfCellZone]
            {
                localX
                {
                    type    ReynoldsPower;
                    coeff   1.84; //
                    exp     -0.2;
                }
                localY
                {
                    type    ReynoldsPower;
                    coeff   1.84; //
                    exp     -0.2;
                }
                localZ
                {
                    type    ReynoldsPower;
                    coeff   0.184; //
                    exp     -0.2;
                }

            }
        }
    }


.. .. toctree::
..    :maxdepth: 1

..    FSHeatTransferCoefficientModels
..    FSdrag




.. note::

    Anisotropic pressure drops can be set using  three different correlations    for the three different local axes. In addition, it is possible to use an
    anisotropic hydraulic diameter ((see :ref:`the structureProperties sub-dictionary <modules_thermalHydraulics_porousMedium_structureProperties>`)). The anisotropy of the hydraulic diameter can
    be set using the keyword *localDhAnisotropy* and assigned to it a vector of
    three scaling factors, one for each local direction.