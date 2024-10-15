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
    Foam::implicitContactFvPatchVectorField

Description
    Frictional contact boundary condition, similar to implicit contact.
    
    The normal contact component of the stress is discretized implicitly. 
    A blending coefficient is used to transition smoothly between the two conditions
    of open and closed gap.

    The tangent traction is always treated explicitly (as in a fixed gradient BC).
    However the details of the implementation change when moving from open gap
    conditiosn (zero tangent traction) to full contact. In particular, in the 
    latter condition, the tangent component of the gradient is added implicitly
    (see gradientBoundaryCoeff and gradientInternalCoeff functions), and subtracted
    explicitly (see explicitTraction in the source code). This helps convergence.
    
    The user can select the amount of friction allowed with the frictionCoefficient.
    If large enough, the material will tend to bond when in contact.
    However, note that the full bonding traction is applied only when the 
    blending coefficient approaches 1 (and the frictionCoefficient is infinite).

    If glued is set to true, the material are assumed to be bonded (blending 
    coefficient is 1).

Usage
    \verbatim
    "fuelOuter|claddingInner"
    {
        type            implicitContact;
        patchType       regionCoupledOFFBEAT;

        // Glue or not the patch
        glued off;
    
        // Contact & dashpot spring moduli
        contactSpringModulus 1e6;
        frictionSpringModulus 1e4;
        dashpotModulus 1;

        frictionCoefficient        1e6;      
        
        relax 0.01;
        
        traction        uniform (0 0 0);
        pressure        uniform 2.25e6;
        value           $internalField;
    }
    \endverbatim
        

SourceFiles
    implicitContactFvPatchVectorField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef implicitContactFvPatchVectorField_H
#define implicitContactFvPatchVectorField_H

#include "regionCoupledBaseOFFBEATFvPatch.H"
#include "LduInterfaceField.H"
#include "fvPatchField.H"
#include "coupledFvPatchField.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class implicitContactFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class implicitContactFvPatchVectorField
:
    public coupledFvPatchField<vector>
{

protected:

    // Private data

        //- Local reference to region couple patch
        const regionCoupledBaseOFFBEATFvPatch& regionCoupledPatch_;

        //- Local patch weights (due to dashpot and local cell)
        mutable scalarField wP_;

        //- Coupled patch weights (due to contact spring, owner and nbr cells)
        mutable scalarField wC_;
        mutable scalarField wC_t_;

        //- Coupled patch weights (due to contact spring, owner and nbr cells)
        mutable scalarField penaltyFact_;
        mutable scalarField penaltyFact_t_;
        
        //- Explicit component of the boudary normal gradient
        mutable vectorField gradExp_;
        
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

        //- Ratio between current friction and stick friction force.
        mutable scalarField frictionToStickRatio_;

        //- Current value of slip between contact patches
        mutable vectorField slip_;

        //- Old value of slip between contact patches
        mutable vectorField slip0_;
        
        //- Friction coefficient
        scalar frictionCoeff_; 

        //- Spring modulus
        scalar contactSpringModulus_; 

        //- Friction spring modulus
        scalar frictionSpringModulus_; 

        //- Dashpot modulus
        scalar dashpotModulus_; 

        //- Takes track of current time value
        scalar currentTime_;
        
        // Relaxation factor for normal traction
        scalar relax_; 

        //- Treat surfaces as attached when they are in contact.
        Switch glued_;  

    // Private functions

        //- Correct weights for this coupled field
        virtual void correctWeights();

        //- Return the gap width
        tmp<scalarField> gapWidth() const;

        //- Return the neighbour patchField
        inline const implicitContactFvPatchVectorField& neighbour() const
        {
            return refCast<const implicitContactFvPatchVectorField>
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


public:

    //- Runtime type information
    TypeName("implicitContact");


    // Constructors

        //- Construct from patch and internal field
        implicitContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        implicitContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given implicitContactFvPatchVectorField
        // onto a new patch
        implicitContactFvPatchVectorField
        (
            const implicitContactFvPatchVectorField&,
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
                new implicitContactFvPatchVectorField(*this)
            );
        }
#endif  

        //- Construct as copy setting internal field reference
        implicitContactFvPatchVectorField
        (
            const implicitContactFvPatchVectorField&,
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
                new implicitContactFvPatchVectorField(*this, iF)
            );
        }


    //- Destructor
    virtual ~implicitContactFvPatchVectorField()
    {}


    // Member functions

        // Access
        
            //- Interface pressure
            inline const scalarField& interfaceP() const
            {
                return interfaceP_;
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
