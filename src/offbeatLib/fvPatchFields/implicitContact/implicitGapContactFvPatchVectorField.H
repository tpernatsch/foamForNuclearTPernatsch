/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::implicitGapContactFvPatchVectorField

Description
    The `implicitGapContact` boundary condition is an extension of the standard 
    `gapContact` BC, introducing an implicit treatment for better stability 
    and convergence. This BC was originally developed to solve complex multi-body 
    contact problems, such as pellet-to-pellet interaction in nuclear fuel rods.

    \note
    A `gapGasModel` is required for this boundary condition to function properly.

    \warning
    This boundary condition assumes that the same boundary condition is applied 
    to a pair of `regionCoupledOFFBEAT` patches, coupled using the AMI (Arbitrary 
    Mesh Interface) mapping algorithm. It is strongly suggested to assign the 
    `owner` role to the smaller surface, e.g., the `fuelOuter` surface in most cases.

    The normal traction is treated implicitly. When no penetration is detected, 
    the boundary condition behaves as a standard `tractionDisplacement` and applies
    the gap pressure. If penetration occurs, a contact pressure proportional to 
    the penetration depth is applied. The penalty factor is calculated similarly
    to the `gapContact` BC, but the resulting traction vector is partly treated
    implicitly.
    
    The tangential traction is treated semi-implicitly up to the maximum friction 
    limit defined by the friction coefficient times the normal pressure. Also the
    friction penalty factor is calculated similarly to the `gapContact` BC.
    The resulting tangent traction vector is partly treated implicitly (scaled by the 
    frictionToStick ratio).

    <b>Augmented Lagrangian Method:</b> 
    To minimize penetration further, an augmented Lagrangian-like approach is 
    available. This method gradually increases the penalty factor over iterations 
    until the penetration falls below a user-defined tolerance. This approach allows 
    smoother convergence, but smaller initial penalty values should be used. 

________________________________________________________________________________

\par Options
    The `implicitGapContact` fvPatchField can be specified in the `boundaryField` 
    subdictionary of the displacement field.

    Parameters in the patch subdictionary for `implicitGapContact`:

    - <b>`patchType`</b> - Specifies the patch field type.  
    **Must be set to** `regionCoupledOFFBEAT`.
    - <b>`penaltyFactor`</b> - Penalty scale factor for normal pressure.  
    **Default:** `1`. **Typical values:** `0.01` to `1`.
    - <b>`penaltyFactorFriction`</b> - Penalty scale factor for frictional forces.  
    **Default:** `1`. **Typical values:** `0.01` to `1`.
    - <b>`frictionCoefficient`</b> - Coefficient of friction between surfaces.  
    **Default:** `0`.
    - <b>`method`</b> - Algorithm for penalty factor adjustment. Either `implicitPenalty`
    (default penalty-based approach) or `augmentedLagrangian`.
    - <b>`relativePenetrationTolerance`</b> - Penetration threshold for the 
    augmented Lagrangian method, defined as a fraction of the cell depth. Set to SMALL
    by default. Te maximum between relative and absolute tolerance is used.
    - <b>`absolutePenetrationTolerance`</b> - Penetration threshold for the 
    augmented Lagrangian method, defined as an absolute value in m. Set to SMALL
    by default. Te maximum between relative and absolute tolerance is used.
    - <b>`relax`</b> - Relaxation factor for implicit contributions.
    **Default:** `1`. **Typical values:** `0.01` to `1`.
    - <b>`glued`</b> - Treats patches as glued. **Default:** `off`.
    - <b>`value`</b> - Initial displacement field value.

________________________________________________________________________________

Usage
    To use the `implicitGapContact` boundary condition, define it in the 
    patch subdictionary inside the `boundaryField` subdictionary.

    Example configuration:

    \verbatim
        "fuelOuter|claddingInner"
        {
            type            implicitGapContact;
            patchType       regionCoupledOFFBEAT;

            glued off;

            // Alternative method
            // method augmentedLagrangian;
            // relativePenetrationTolerance 0.001;
            // absolutePenetrationTolerance 1e-7;

            penaltyFactor           0.1;
            penaltyFactorFriction   0.001;
            penaltyFactorDashpot    0.0;

            frictionCoefficient     0.1;

            relax 0.1;

            value           $internalField;
        }
    \endverbatim        

SourceFiles
    implicitGapContactFvPatchVectorField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef implicitGapContactFvPatchVectorField_H
#define implicitGapContactFvPatchVectorField_H

#include "regionCoupledBaseOFFBEATFvPatch.H"
#include "LduInterfaceField.H"
#include "fvPatchField.H"
#include "coupledFvPatchField.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class implicitGapContactFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class implicitGapContactFvPatchVectorField
:
    public coupledFvPatchField<vector>
{

protected:

    // Private data

        //- Local reference to region couple patch
        const regionCoupledBaseOFFBEATFvPatch& regionCoupledPatch_;

        //- Coupled patch weights (due to contact spring, owner and nbr cells)
        mutable scalarField wC_;
        mutable scalarField wC_t_;

        //- Coupled patch weights (due to contact spring, owner and nbr cells)
        mutable scalarField penaltyFact_;
        mutable scalarField penaltyFact_t_;
        
        //- Explicit component of the boudary displacement
        mutable vectorField valueExp_;
        
        //- Explicit component of laplacian
        mutable vectorField laplacianExp_; 
        
        // Surface pressure
        scalarField pressure_;
        
        // Interfacial Pressure
        mutable scalarField interfaceP_;
        
        //- Gap widht
        mutable scalarField gapWidth_;
        
        //- Keep track if gapWidth has been updated in this or nbr patch
        //  in current iteration
        mutable bool gapCalculated_;

        //- Ratio between current friction and stick friction force.
        mutable scalarField frictionToStickRatio_;

        //- Current value of slip between contact patches
        mutable vectorField slip_;

        //- Old value of slip between contact patches
        mutable vectorField slip0_;
        
        //- Friction coefficient
        scalar frictionCoeff_; 

        //- Contact penalty scale factor
        scalar penaltyScaleFact_; 

        //- Contact penalty scale factor for friction
        scalar penaltyScaleFactFriction_; 

        //- Contact penalty scale factor for dashpot
        scalar penaltyScaleFactDashpot_; 

        //- Takes track of current time value
        scalar currentTime_;
        
        // Relaxation factor for normal traction
        scalar relax_; 

        //- Treat surfaces as attached when they are in contact.
        Switch glued_;  

        //- Calculate implicit component using an implicit penalty approach
        // or an augmentedLagrangian approach
        word method_;

        //- Relative penetration tolerance for augmentedLagrangian approach
        scalar relTolerance_; 

        //- Absolute penetration tolerance for augmentedLagrangian approach
        scalar absTolerance_; 

        //- Offset for gapWidth (a positive value increases the gap width)
        scalar offset_;

        //- Indicates whether patch field is coupled (during construction and solution)
        // or not coupled (after solve)
        mutable bool coupled_;

    // Private functions

        //- Correct weights for this coupled field
        virtual void correctWeights();

        //- Return the gap width
        tmp<scalarField> gapWidth() const;

        //- Return the neighbour patchField
        inline const implicitGapContactFvPatchVectorField& neighbour() const
        {
            return refCast<const implicitGapContactFvPatchVectorField>
            (
                regionCoupledPatch_.neighbFvPatch().lookupPatchField
                <volVectorField, vector>
                (
                    internalField().name()
                )
            );
        }

        //- Minimum patch stiffness for interface pressure evaluation
        scalarField boundaryStiffness() const;

        //- Minimum patch shear stiffness for friction evaluation
        scalarField boundaryShearStiffness() const;


public:

    //- Runtime type information
    TypeName("implicitGapContact");


    // Constructors

        //- Construct from patch and internal field
        implicitGapContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        implicitGapContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given implicitGapContactFvPatchVectorField
        // onto a new patch
        implicitGapContactFvPatchVectorField
        (
            const implicitGapContactFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

#ifdef OPENFOAMESI
         //- Construct and return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return tmp<fvPatchField<vector>>
            (
                new implicitGapContactFvPatchVectorField(*this)
            );
        }
#endif  

        //- Construct as copy setting internal field reference
        implicitGapContactFvPatchVectorField
        (
            const implicitGapContactFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new implicitGapContactFvPatchVectorField(*this, iF)
            );
        }


    //- Destructor
    virtual ~implicitGapContactFvPatchVectorField()
    {}


    // Member functions

        // Access
        
            //- Interface pressure
            inline const scalarField& interfaceP() const
            {
                return interfaceP_;
            } 

            //- Return true if this patch field is derived from
            //  coupledFvPatchField<Type>.
            virtual bool coupled() const
            {
                return coupled_;
            }

            //- Make field coupled or uncoupled
            virtual void coupleField(bool coupleOrUnCouple)
            {
                coupled_ = coupleOrUnCouple;
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

            //- Return neighbour coupled internal cell data
            virtual tmp<vectorField> patchNeighbourField() const;

            //- Retrieve neighbour coupled data
            virtual void patchNeighbourField(UList<vector>&) const;

            //- Update the coefficients associated with the patch field
            //  Sets Updated to true
            virtual void updateCoeffs();
            
            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType
            );

            //- Return patch-normal gradient
            virtual tmp<vectorField> snGrad() const;

            //- Return patch-normal gradient
            //  Note: the deltaCoeffs supplied are not used
            virtual tmp<vectorField> snGrad
            (
                const scalarField& deltaCoeffs
            ) const;
            
            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<vectorField> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<vectorField> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<vectorField> gradientInternalCoeffs
            (
                const scalarField& deltaCoeffs
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<vectorField> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<vectorField> gradientBoundaryCoeffs
            (
                const scalarField& deltaCoeffs
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<vectorField> gradientBoundaryCoeffs() const;
            
        // Coupled interface functionality

#ifdef OPENFOAMFOUNDATION
            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<scalar>& result,
                const scalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<vector>&,
                const Field<vector>&,
                const scalarField&,
                const Pstream::commsTypes commsType
            ) const;
#elif OPENFOAMESI
             //- Update result field based on interface functionality
             virtual void updateInterfaceMatrix
             (
                 solveScalarField& result,
                 const bool add,
                 const lduAddressing& lduAddr,
                 const label patchId,
                 const solveScalarField& psiInternal,
                 const scalarField& coeffs,
                 const direction,
                 const Pstream::commsTypes commsType
             ) const;
  
             //- Update result field based on interface functionality
             virtual void updateInterfaceMatrix
             (
                 Field<vector>&,
                 const bool add,
                 const lduAddressing& lduAddr,
                 const label patchId,
                 const Field<vector>&,
                 const scalarField&,
                 const Pstream::commsTypes commsType
             ) const;            
#endif              

            virtual void manipulateMatrix(fvMatrix<vector>& matrix);      


            //- Return face-gradient transform diagonal
            virtual tmp<Field<vector>> snGradTransformDiag() const;

            //- Return face-gradient tangential transform diagonal
            virtual tmp<Field<vector>> tGradTransformDiag() const;


            //- Return the interface type
            virtual const word& interfaceFieldType() const
            {
                return regionCoupledPatch_.regionCoupleType();
            }


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
