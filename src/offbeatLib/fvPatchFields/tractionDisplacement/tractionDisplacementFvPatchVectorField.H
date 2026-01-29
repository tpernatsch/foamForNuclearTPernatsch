/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::tractionDisplacementFvPatchVectorField

Description
    The `tractionDisplacement` boundary condition is used to apply a fixed or time-
    dependent traction vector. The total traction vector is composed 
    of a user-defined `traction` (provided as a vector) and a `pressure`, assumed 
    positive when acting against the surface normal vector.

    Optionally, a plane strain approximation can be activated: it imposes the same 
    uniform normal strain on the patch. It can be useful when simulating 2D discrete 
    rods with coarse axial division.
    
    Finally, a spring-dashpot system can be activated: it simulates the presence of 
    a spring and dashpot attached to the patch, and it can help convergence 
    for bodies that are not well constrained.
    
    The spring-dahpot system produce a force

    $$
   \begin{aligned}
   \tau_{sd} \cdot A_f = -K_s \cdot (D_f) -K_d \cdot (D_f - D_{f,prev})
   \end{aligned}
   $$

    where:

    - \f$\tau_{sd}\f$ is the spring-dashpot force on a patch face,
    - \f$A_f\f$ is the face area,
    - \f$K_s\f$ is the spring modulus,
    - \f$K_d\f$ is the dashpot modulus,
    - \f$D_f\f$ is the current face displacement vector, and
    - \f$D_{f, prev}\f$ is the face displacement vector at the previous iteration.
    
    From this expression, it follows that the dashpot force vanishes at convergence, 
    meaning that the dashpot acts only during the iterative process to stabilize 
    and aid in achieving a solution. The spring force, however, persists and contributes 
    to the equilibrium state.

________________________________________________________________________________

\par Options
    The `tractionDisplacement` fvPatchField can be selected in the patch subdictionary
    inside the `boundaryField` subdictionary of the displacement field.

    Parameters in the patch subdictionary for `tractionDisplacement`:

    - <b>`tractionList`</b> - A time-dependent traction table. Each entry consists of 
    a time step and a corresponding traction vector value in Pa.
    - <b>`traction`</b> - Fixed traction vector in Pa. Used only if `tractionList` 
    is absent.
    - <b>`pressureList`</b> - A time-dependent pressure table. Each entry consists of 
    a time step and a corresponding pressure value in Pa.
    - <b>`pressure`</b> - Fixed pressure in Pa. Used only if `pressureList` is absent.
    - <b>`planeStrain`</b> - Activates the plane strain approximation for the normal 
    stress at the boundary. When enabled, the normal component of the strain (e.g. epislonZZ)
    is assumed constant at the boundary. Useful for long axisymmetric rods. It cannot
    be used in combination with `flatSurface` keyword. **Default: `false`.**
    - <b>`flatSurface`</b> - Activates the flat surface approximation for the normal 
    stress at the boundary. When enabled, the normal component of the displacement (e.g. Dz)
    is assumed constant at the boundary, i.e. the patch is forced to remain flat. 
    Useful for instance for the simulation of a half pellet fragment. It cannot
    be used in combination with `planeStrain` keyword. **Default: `false`.**
    - <b>`fixedSpring`</b> - Activates a fixed spring-dashpot system for additional 
    stability. **Default: `false`.**
    - <b>`fixedSpringModulus`</b> - Spring modulus in N/m. Required when `fixedSpring` 
    is set to `true`.
    - <b>`dashpotModulus`</b> - Dashpot modulus in N/m. Required when `fixedSpring` 
    is set to `true`.
    - <b>`relax`</b> - Relaxation factor for gradient updates. **Default: `1.0`.**
    - <b>`value`</b> - Initial displacement value (not stress).

________________________________________________________________________________

Usage
    To apply the `tractionDisplacement` boundary condition to a given patch, the 
    following examples can be used as templates:

    Example 1: Time-varying traction and pressure

    \verbatim
        cladOuter
        {
            type            tractionDisplacement;

            tractionList
            {
                type    table;

                values
                (
                    (0.0 (0 0 0))
                    (100.0 (-1000 0 0))
                );

                outOfBounds         clamp;     // Optional out-of-bounds handling
                interpolationScheme linear;    // Optional interpolation method
            }

            pressureList
            {
                type    table;

                values
                (
                    (0.0    1e5)
                    (100.0  1e7)
                );

                outOfBounds         clamp;
                interpolationScheme linear;
            }

            planeStrain     false;  // Optional: Plane strain approximation
            flatSurface     false;  // Optional: Flat surface approximation
            relax           1.0;    // Relaxation factor
            value           $internalField;
        }
    \endverbatim

    Example 2: Constant pressure and zero traction

    \verbatim
        cladOuter
        {
            type            tractionDisplacement;

            traction        uniform (0 0 0);  // Fixed traction
            pressure        uniform 15.5e6;  // Fixed pressure

            planeStrain     false;  // Optional: Plane strain approximation
            flatSurface     false;  // Optional: Flat surface approximation

            relax           1.0;
            value           $internalField;
        }
    \endverbatim

    Example 3: Adding a fixed spring-dashpot system

    \verbatim
        cladOuter
        {
            type            tractionDisplacement;

            traction        uniform (0 0 0);  // Fixed traction
            pressure        uniform 15.5e6;  // Fixed pressure
            
            planeStrain     false;  // Optional: Plane strain approximation
            flatSurface     false;  // Optional: Flat surface approximation
            
            relax           1.0;

            fixedSpring     on;              // Activate spring-dashpot system
            fixedSpringModulus 1e2;          // Spring modulus in N/m
            dashpotModulus     1e4;          // Dashpot modulus in N/m

            value           $internalField;
        }
    \endverbatim

    \note
    You can mix time-dependent and fixed values for traction and pressure in the 
    boundary condition definition.

    \warning
    For incremental simulations (`DD` instead of `D`), the required traction and 
    pressure values in the dictionary are the respective total quantities (not 
    incremental), even for time-dependent cases.

SourceFiles
    tractionDisplacementFvPatchVectorField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef tractionDisplacementFvPatchVectorField_H
#define tractionDisplacementFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidDirectionMixedFvPatchVectorField.H"
#include "Table.H"

#ifdef OPENFOAMESI
    #include "Function1.H"
#endif    


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class tractionDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class tractionDisplacementFvPatchVectorField
:
    public solidDirectionMixedFvPatchVectorField
{

    // Private Data
        
protected:

    // Protected Data
  
        //- Surface traction vector
        vectorField traction_;
        
        //- Surface pressure
        scalarField pressure_;
  
        //- Total traction vector (pressure + traction)
        vectorField totalTraction_;
        
        //- Old total traction vector (pressure + traction)
        //  Only for incremental solvers
        autoPtr<vectorField> totalTraction0_;
        
        //- relaxation factor for update
        scalar relax_;

        //- Name of stress field
        word stressName_;   

        //- Apply plane strain correction 
        //  (normal component of normal gradient kept uniform)
        Switch planeStrain_;

        //- Apply flat surface correction 
        //  (normal component of displacement kept uniform)
        Switch flatSurface_;

#ifdef OPENFOAMFOUNDATION
        //- Time varying list of pressure values
        autoPtr<Function1s::Table<scalar> > pressureList_;

        //- Time varying list of traction values
        autoPtr<Function1s::Table<vector> > tractionList_;

#elif OPENFOAMESI
        typedef Foam::Function1Types::Table<scalar> scalarTable;
        typedef Foam::Function1Types::Table<vector> vectorTable;

        //- Time varying list of pressure values
        autoPtr<scalarTable > pressureList_;

        //- Time varying list of traction values
        autoPtr<vectorTable >  tractionList_;
#endif        

        //- Apply fixed spring-dashpot system or not
        Switch fixedSpring_;

        //- Spring modulus
        scalar fixedSpringModulus_; 

        //- Dashpot modulus
        scalar dashpotModulus_; 

        // //- Patch weights for spring-dashpot system
        // scalarField w_;

        //- Value of current time step (in s)
        mutable scalar currentTime_;
        
        //- Total force on patch
        vector Ftot_;
        
public:

    //- Runtime type information
    TypeName("tractionDisplacement");


    // Constructors

        //- Construct from patch and internal field
        tractionDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary.
        //  The bool is used to decide whether to initialize the fields
        tractionDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&,
            const bool valueRequired = true
        );

        //- Construct by mapping given
        //  tractionDisplacementFvPatchVectorField onto a new patch
        tractionDisplacementFvPatchVectorField
        (
            const tractionDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );

        //- Construct as copy setting internal field reference
        tractionDisplacementFvPatchVectorField
        (
            const tractionDisplacementFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new tractionDisplacementFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            virtual const vectorField& traction() const
            {
                return traction_;
            }

            virtual vectorField& traction()
            {
                return traction_;
            }

            virtual const scalarField& pressure() const
            {
                return pressure_;
            }

            virtual  scalarField& pressure()
            {
                return pressure_;
            }

            inline const vector& totalForce() const
            {
                return Ftot_;
            }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        // Evaluation functions

            // //- Return gradient at boundary
            // virtual tmp<Field<vector>> snGrad() const;

            //- Update the pressure and traction
            virtual void updateTraction();

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            // //- Evaluate the patch
            // virtual void evaluate(const Pstream::commsTypes);
 
            // //- Return the matrix diagonal coefficients corresponding to the
            // //  evaluation of the gradient of this patchField
            // virtual tmp<Field<vector>> gradientInternalCoeffs() const;
 
            // //- Return the matrix source coefficients corresponding to the
            // //  evaluation of the gradient of this patchField
            // virtual tmp<Field<vector>> gradientBoundaryCoeffs() const; 
 
            // //- Return the matrix diagonal coefficients corresponding to the
            // //  evaluation of the value of this patchField with given weights
            // virtual tmp<Field<vector>> valueInternalCoeffs
            // (
            //     const tmp<scalarField>&
            // ) const;

            // //- Return the matrix source coefficients corresponding to the
            // //  evaluation of the value of this patchField with given weights
            // virtual tmp<Field<vector>> valueBoundaryCoeffs
            // (
            //     const tmp<scalarField>&
            // ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
