/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::rheologyByMaterial

Description
    Rheology or constitutive laws are specified on each material subdictionary
    within the `solverDict`. This class also includes some functions that allow 
    OFFBEAT to run simplified simulations employing the plane-stress (e.g. for 
    fuel r-theta dics) or the standard 1.5D approximation.

________________________________________________________________________________

\par Options
    A few options can be selected in the `rheologyOptions` subdictionary located
    within the `solverDict`. 

    Parameters in `rheologyOptions`:

    - <b>`planeStress`</b> - If activated, the strain tensor z-component is 
    modified to consider the plane stress approximation. **It is off by default**. 
    - <b>`modifiedPlaneStrain`</b> - If activated, the strain tensor z-component 
    is modified to consider the typical modified-plane strain approximation of standard 
    fuel performance. **It is off by default**. 
    - <b>`precisionSpring`</b> - Only valid if `modifiedPlaneStrain` is activated, 
    it sets the convergence criterion for the modified plane strain loop. **Set to 
    1e-6 by default**.
    - <b>`springModulus`</b> - Only valid if `modifiedPlaneStrain` is activated, 
    spring modulus in N/m for the plenum spring. **Set to 1e-6 by default**. 
    - <b>`springPreCompression`</b> - Only valid if `modifiedPlaneStrain` is activated, 
    pre-compression in m for the plenum spring. **Set to 0 by default**. 
    - <b>`coolantPressureList`</b> - The coolant pressure that acts on the top 
    cladding cap. 
    - <b>`thermalExpansion`</b> - Activate or deactivate thermal strains. 
    **It is on by default**.
    - <b>`solvePressureEqn`</b> - If activated, an additional equation is 
    solved for the hydrostatic stress (`sigmaHyd`) in order to quell or smooth 
    some oscillations that can appear in the stress field, in particular in cases 
    of considerable plastic strains. **It is off by default**
    - <b>`pressureSmoothingFactor`</b> - Only read if `solvePressureEqn` is 
    activated, it sets the smoothing scale factor for the `sigmaHyd` equation. 
    **It is 1.0 by default**.

    \warning
    When activating the `solvePressureEqn` keyword an additional equation is solved 
    for the field `sigmaHyd`. This means that a linear solver for this field should 
    be added to the `fvSolution` file in the `system` folder (as it is typically done 
    for `D` or `T`).
    <br><br>
    Similarly a `sigmaHyd` file must be present in the initial folder (e.g. `0/` folder). 
    The boundary condition for all non-empty or non-wedge patches should be of the type 
    `extrapolatedFixedValue`.
    <br><br>
    Finally, for standard fuel performance simulations it is not suggested to use the 
    pressure smoothing equation option. If any oscillation arises in the stress field, 
    it should be solved either by refining the mesh, by selecting appropriate relaxation
    factors or by using the `multiMaterialCorrection` in the settings of `mechanicsSolver`.
    
    \warning
    The `planeStress` and `modifiedPlaneStrain` options are valid only for models 
    with empty boundary conditions along the z-direction.
    <br><br>
    Additionally, `planeStress` has priority over `modifiedPlaneStrain` i.e. 
    `modifiedPlaneStrain` can be activated only if  `planeStress` switch is not activated.

________________________________________________________________________________

Usage
    Here is a snippet of the `solverDict` that shows how to correctly set rheology 
    class for a standard simulation (without planeStrain or planeStress approximation):

    \verbatim
        rheology byMaterial;

        rheologyOptions
        {
            // thermalExapansion is on by default
            thermalExpansion    on;

            // plane stress and plane strain are off by default
            planeStress         off;
            modifiedPlaneStrain off;

            // smoothing of pressure equation is off by default
            solvePressureEqn off;
        }

        // List of materials, one per cellZone.
        materials
        {
            fuel
            {
                material UO2;
                ...

                rheologyModel elasticity;
            }

            cladding
            {
                material zircaloy;
                ...

                rheologyModel misesPlasticity;

                ...
            }

            // ... sub-dictionaries for remaining materials, if any
        }
    \endverbatim

    Here instead is an example showing how one can set the `rheologyOptions` 
    sub-dictionary inside the `solverDict` file to activate the modifiedPlaneStrain 
    approximation (i.e. **for a 1.5D simulation**):

    \verbatim
        rheologyOptions
        {
            // thermalExapansion is on by default
            thermalExpansion    on;

            // plane stress and plane strain approximation are off by default
            planeStress         off;
            modifiedPlaneStrain on;

            // In case of plane strain, the user might specify the following:
            // - convergence criterion for spring force balance  
            // - spring modulus
            // - pre-compression (initial compression, taken as positive in [m])
            // - the coolant pressure that acts on the top cladding cap
            precisionSpring      1e-6;
            springModulus        3.5e3;       
            springPreCompression 0.0;

            coolantPressureList
            {
                type                table;

                // Instead of "file", the user can insert "values"
                file                "$FOAM_CASE/constant/lists/topCladdingRingPressureCoolantList";

                // values
                // (
                //     (0.0 1e5)
                //     (100.0 1e6)
                // );

                format              foam;      // data format (optional)
                outOfBounds         clamp;     // optional out-of-bounds handling
                interpolationScheme linear;    // optional interpolation method
            }

            // smoothing of pressure equation is off by default
            solvePressureEqn off;
        }
    \endverbatim

SourceFiles
    rheologyByMaterial.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef rheologyByMaterial_H
#define rheologyByMaterial_H

#include "rheology.H"
#include "sliceMapper.H"
#include "Table.H"

#ifdef OPENFOAMESI
    #include "Function1.H"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class rheologyByMaterial Declaration
\*---------------------------------------------------------------------------*/

class rheologyByMaterial
:
    public rheology
{
    // Private data
        
        //- Precision threshold for spring+friction force balance 
        scalar precision_;

        //- Spring modulus
        scalar springModulus_;

        //- Spring initial loading (in m)
        scalar springPreCompression_;

        //- Pointer to slice mapper (necessary for modifiedPlaneStrain)
        const sliceMapper* mapper_;

        //- Friction coefficient
        scalar frictionCoeff_;

        //- Contact penalty scale factor for friction
        scalar penaltyScaleFactFriction_; 

        //- Normal pressure below which friction is ramped linearly; 
        //  above it, full Coulomb friction is applied.
        scalar frictionSaturationPressure_;

        //- Time tracker
        scalar timeCheck_;  

        //- Relaxation factor for spring+friction force in modifiedPlaneStrain
        scalar relax_;     

#ifdef OPENFOAMFOUNDATION
        //- Fluid pressure time-dependent list (necessary for modifiedPlaneStrain)
        autoPtr<Function1s::Table<scalar>> fluidPressureList_;
#elif OPENFOAMESI
        //- Fluid pressure time-dependent list (necessary for modifiedPlaneStrain)
        autoPtr<Function1<scalar>> fluidPressureList_;
#endif

        //- SpringForce (only for modifiedPlaneStrain)
        autoPtr<scalarField> springForce_;

        //- SpringForce (only for modifiedPlaneStrain)
        autoPtr<scalarField> frictionForce_;

        //- Pressure force (only for modifiedPlaneStrain)
        autoPtr<scalarField> pressureForce_;

        //- Old time slip (only for modifiedPlaneStrain)
        autoPtr<scalarField> slip0_;

        //- Current time slip (only for modifiedPlaneStrain)
        autoPtr<scalarField> slip_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        rheologyByMaterial(const rheologyByMaterial&);

        //- Disallow default bitwise assignment
        void operator=(const rheologyByMaterial&);

public:

    //- Runtime type information
    TypeName("byMaterial");
    

    // Constructors

        //- Construct from mesh, materials and dictionary
        rheologyByMaterial
        (
            const fvMesh& mesh,
            materials& mat,
            const dictionary& mechanicsDict
        ); 
    

    //- Destructor
    ~rheologyByMaterial();


    // Member Functions

        //- Correct/update the properties
        virtual void correct
        (
            volSymmTensorField& sigma, 
            volSymmTensorField& epsilon,
            volVectorField& D
        );     

        //- Correct out of plane strain for planeStress cases
        virtual void correctPlaneStress
        (
            volSymmTensorField& epsilon
        );      

        //- Correct out of plane strain for modifiedPlaneStress cases
        virtual void correctModifiedPlaneStrain
        (
            volSymmTensorField& epsilon,
            volVectorField& D
        );     

        //- Update total fields at the end of time step
        virtual void updateTotalFields();

        //- Collect additional strain field from zone constitutive laws
        volSymmTensorField getLawsAdditionalStrain();

        //- Return next deltaT according to rheology
        virtual scalar nextDeltaT(); 

        //- Return next macroTime according to rheology
        virtual scalar nextMacroTime(); 
          
        //- Return maximum time according to rheology
        virtual scalar lastMacroTime();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
