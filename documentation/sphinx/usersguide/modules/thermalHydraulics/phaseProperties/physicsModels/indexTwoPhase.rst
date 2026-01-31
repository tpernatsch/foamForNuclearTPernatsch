.. _modules_thermalHydraulics_porousMedium_physicsModelsTwoPhase:

The *physicsModels* sub-dictionary for two-phase flows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The *physicsModels* sub-dictionary is used in two-phase simulations to define the models that describe the interaction between fluid1 and structure, fluid2 and structure, as well as between fluid1 and fluid 2. Several types of models can be defined:

- :ref:`fluid-fluid drag models <modules_thermalHydraulics_porousMdeium_FFdrag>`  and :ref:`fluid-structure drag models <modules_thermalHydraulics_porousMdeium_FSdrag>`  via the ``dragModels`` sub-dictionary;
- :ref:`two-phase drag multiplier models <modules_thermalHydraulics_porousMdeium_multipliers>`  via the ``twoPhaseDragMultiplierModel`` sub-dictionary;
- :ref:`virtual mass coefficient models <modules_thermalHydraulics_porousMedium_virtualMassCoefficientModels>` via the ``virtualMassCoefficientModel`` sub-dictionary;
- pair geometry models via the ``pairGeometryModels`` sub-dictionary, which in turn includes;
   - :ref:`dispersion models <modules_thermalHydraulics_porousMedium_dispersionModels>` for fluid-fluid pairs via the ``dispersionModel`` sub-dictionary;
   - :ref:`interfacial area density models <modules_thermalHydraulics_porousMedium_interfacialAreaModels>` for fluid-fluid pairs via the ``interfacialAreaDensityModel`` sub-dictionary;
   - :ref:`contact partition models <modules_thermalHydraulics_porousMedium_contactPartitionModels>` for fluid-structure pairs via the ``contactPartitionModel`` sub-dictionary
- :ref:`phase change models <modules_thermalHydraulics_porousMedium_phaseChangeModels>` via the ``phaseChangeModel`` sub-dictionary, which in turn includes:
   - :ref:`latent heat models <modules_thermalHydraulics_porousMedium_latentHeatModels>` for fluid-fluid pairs via the ``latentHeatModel`` sub-dictionary;
   - :ref:`saturation models <modules_thermalHydraulics_porousMedium_saturationModels>` for fluid-structure pairs via the ``saturationModel`` sub-dictionary


Each of these two sub-dictionaries consists of a series of other sub-dictionaries whose properties are applied to the cell zones that have their name in the sub-dictionary keys.  Three options exist :
- Use ``constant`` as  type and give a constant value using the keyword ``value``.
- Use as type a single model chosen among those available;
- Use ``byRegime`` as  type, and then chose one of the available models for each regime. This requires having set a :ref:`regime map <modules_thermalHydraulics_porousMedium_regimeMapModels>`.


An example of dictionary is reported below. 

.. code :: cpp

   physicsModels
   {
      dragModels
      {
         "liquid.vapour"
         {
               type    SchillerNaumann;
         }

         "liquid.structure"
         {
               "low:mid:top"
               {
                  type    ReynoldsPower;
                  coeff   2; // coeff A in A*Re^b
                  exp     -0.125; // coeff b in A*Re^b
               }
         }
      }

      twoPhaseDragMultiplierModel
      {
         type                Kaiser88;
         multiplierFluid     liquid;
      }

      heatTransferModels
      {
         "liquid.vapour"
         {
               "liquid"
               {
                  type    NusseltReynoldsPrandtlPower;
                  const   10;
                  coeff   0;
                  expRe   0;
                  expPr   0;
               }
               "vapour"
               {
                  type    NusseltReynoldsPrandtlPower;
                  const   10;
                  coeff   0;
                  expRe   0;
                  expPr   0;
               }
         }

         "liquid.structure"
         {
               "mid"
               {
                  type                superpositionNucleateBoiling;
                  forcedConvection
                  {
                     type            NusseltReynoldsPrandtlPower;
                     const           7.48467;
                     coeff           0.02994;
                     expRe           0.77;
                     expPr           0.77;
                  }
                  poolBoiling
                  {
                     type                    Shah;
                     useExplicitHeatFlux     false;
                  }
                  flowEnhancementFactor
                  {
                     type            COBRA-TF;
                  }
                  suppressionFactor
                  {
                     type            COBRA-TF;
                  }
               }
         }
      }

      virtualMassCoefficientModel
      {
         type    constant;
         value   0.1;
      }

      pairGeometryModels
      {
         "liquid.vapour"
         {
               dispersionModel
               {
                  type        byRegime;
                  regimeMap   "slugMist";
                  "slug"
                  {
                     type            constant;
                     dispersedPhase  "vapour";
                  }
                  "mist"
                  {
                     type            constant;
                     dispersedPhase  "liquid";
                  }
               }

               interfacialAreaDensityModel
               {
                  type        spherical;
               }
         }

         "liquid.structure"
         {
               contactPartitionModel
               {
                  type        byRegime;
                  regimeMap   "slugMist";
                  "slug"
                  {
                     type        constant;
                     value       1.0;
                  }
                  "mist"
                  {
                     type        constant;
                     value       0.1;
                  }
               }
         }
      }

      phaseChangeModel
      {
         type                heatDriven;
         mode                conductionLimited;
         correctLatentHeat   false;

         latentHeatModel
         {
               type        FinkLeibowitz;
               adjust      true;
         }

         saturationModel
         {
               type        BrowningPotter;
         }
      }
   }




.. The main run-time selctable models include:

.. .. toctree::
..    :maxdepth: 1

..    FFHeatTransferCoefficientModels
..    FSHeatTransferCoefficientModels
..    FFdrag
..    FSdrag
..    multipliersDrag
..    phaseChangeModels
..    contactPartitionModels
..    dispersionModels
..    fluidDiameterModels
..    interfacialAreaModels
..    virtualMassCoefficientModels
..    latentHeatModels
..    saturationModels



.. note::

    Anisotropic pressure drops can be set by by setting three different correlations    for the three different local axes (see :ref:`the physicsModels sub-dictionary <modules_thermalHydraulics_porousMedium_physicsModels>`). In addition, it is possible to use an
    anisotropic hydraulic diameter ((see :ref:`the structureProperties sub-dictionary <modules_thermalHydraulics_porousMedium_structureProperties>`)). The anisotropy of the hydraulic diameter can
    be set using the keyword *localDhAnisotropy* and assigned to it a vector of
    three scaling factors, one for each local direction.