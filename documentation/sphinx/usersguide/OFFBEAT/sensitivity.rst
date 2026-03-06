Sensitivity analysis and uncertainty quantification
===================================================

OFFBEAT provides built-in capabilities for sensitivity analysis (SA) and
uncertainty quantification (UQ) by allowing users to systematically perturb
input parameters.

This is handled through the dictionary:

* ``system/userParameters``

which can be used to modify material properties and behavioural model
parameters.

-------------------------------------------------------------------------------

Input parameter perturbation
----------------------------

The ``userParameters`` dictionary defines:

* scale factors ``F``
* offset values ``delta``

used to perturb the nominal value :math:`v` computed by specific material
properties and behavioural models.

The perturbed value :math:`v'` is computed as:

.. math::

   v' = v \cdot F + \delta

Perturbation factors are defined in ``system/userParameters`` according to the
following structure.

Material properties
^^^^^^^^^^^^^^^^^^^

Perturbations defined in the ``materials`` dictionary refer to thermo-mechanical
properties (e.g. density, emissivity, Young's modulus, etc.).

By default, perturbations set at this level apply to all materials.

To apply perturbations to a specific material (e.g. ``UO2``), define overrides
in a subdictionary such as ``materials/UO2``.

Behavioural models
^^^^^^^^^^^^^^^^^^

Perturbations defined in the ``behavioralModels`` dictionary refer to
behavioural models (e.g. relocation, densification, swelling, etc.).

Correlation-specific perturbations can be defined using an additional
subdictionary. For example, to perturb only the ``UO2Frapcon`` correlation
for swelling, define the entry under:

* ``behavioralModels/swelling/UO2Frapcon``

.. note::

   A comprehensive list of all available perturbation parameters is printed
   in the OFFBEAT log at the beginning of the simulation.

-------------------------------------------------------------------------------

Example: ``system/userParameters``
----------------------------------

.. code-block:: cpp

   /*--------------------------------*- C++ -*----------------------------------*\
     =========                 |
     \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
      \\    /   O peration     | Website:  https://openfoam.org
       \\  /    A nd           | Version:  9
        \\/     M anipulation  |
   \*---------------------------------------------------------------------------*/
   FoamFile
   {
       version         2;
       format          ascii;
       class           dictionary;
       location        "system";
       object          userParameters;
   }

   materials
   {
       // Perturbations defined here impact each property model across all
       // materials defined in the simulation
       F_rho           [ 0 0 0 0 0 ] 1;
       F_Cp            [ 0 0 0 0 0 ] 1;
       F_k             [ 0 0 0 0 0 ] 1;
       F_emissivity    [ 0 0 0 0 0 ] 1;
       F_E             [ 0 0 0 0 0 ] 1;
       F_nu            [ 0 0 0 0 0 ] 1;
       F_alphaT        [ 0 0 0 0 0 ] 1;

       // UO2-specific perturbations
       UO2
       {
           F_rho           [ 0 0 0 0 0 0 0 ] 1.05704;
           F_Cp            [ 0 0 0 0 0 0 0 ] 1.04758;
           F_k             [ 0 0 0 0 0 0 0 ] 0.973783;
           F_emissivity    [ 0 0 0 0 0 0 0 ] 0.904494;
           F_E             [ 0 0 0 0 0 0 0 ] 1.02771;
           F_nu            [ 0 0 0 0 0 0 0 ] 0.951355;
           F_alphaT        [ 0 0 0 0 0 0 0 ] 1.00606;
       }

       // Zircaloy-specific perturbations
       Zircaloy
       {
           F_rho           [ 0 0 0 0 0 0 0 ] 1.03358;
           F_Cp            [ 0 0 0 0 0 0 0 ] 1.08002;
           F_k             [ 0 0 0 0 0 0 0 ] 0.960583;
           F_emissivity    [ 0 0 0 0 0 0 0 ] 0.974784;
           F_E             [ 0 0 0 0 0 0 0 ] 1.05035;
           F_nu            [ 0 0 0 0 0 0 0 ] 0.952649;
           F_alphaT        [ 0 0 0 0 0 0 0 ] 0.905154;
       }
   }

   behavioralModels
   {
       densification
       {
           F_epsilonDensification [ 0 0 0 0 0 0 0 ] 1.08;
       }

       relocation
       {
           F_epsilonRelocation [ 0 0 0 0 0 0 0 ] 0.98;
       }

       swelling
       {
           // Correlation-specific perturbation
           UO2Frapcon
           {
               F_epsilonSwelling [ 0 0 0 0 0 0 0 ] 1.25;
           }
       }
   }