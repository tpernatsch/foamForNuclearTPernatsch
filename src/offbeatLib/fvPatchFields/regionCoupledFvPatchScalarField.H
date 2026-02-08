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
    Foam::regionCoupledFvPatchScalarField

Description
    General region-coupled patchField.
    Conceptually similar to cyclicAMI, but does not assume continuity
    across the interface, i.e. the fvPatch is treated as uncoupled
    from the delta point of view. This allows for exchange of different
    physics across the interface.
    If the patchField is defined as coupled, the fvPatch is an interface
    and is incorporated into the matrix implicitly.
    If not coupled, then the patchField reverts to a simple
    zeroGradient.

Usage
    \verbatim
    fuelOuter
    {
        type            regionCoupled;
        patchType       regionCoupledOFFBEAT;
        coupled         true;
        value           $internalField;
    } 
    \endverbatim

SourceFiles
    regionCoupledFvPatchScalarField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledFvPatchScalarField_H
#define regionCoupledFvPatchScalarField_H

#include "regionCoupledBaseOFFBEATFvPatch.H"
#include "LduInterfaceField.H"
#include "fvPatchField.H"
#include "coupledFvPatchField.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class regionCoupledFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledFvPatchScalarField
:
    public coupledFvPatchField<scalar>
{

protected:

    // Private data

        //- Local reference to region couple patch
        const regionCoupledBaseOFFBEATFvPatch& regionCoupledPatch_;

        //- True if patch is coupled. False reverts to zeroGradient field.
        bool coupled_;
        
        //- Coupled patch weights
        scalarField weights_;
        
        //- Area-weighting correction
        scalarField correctionFactor_;

    // Private functions

        //- Local weight for this coupled field
        virtual tmp<scalarField> weights() const;

        //- Local weight for this coupled field
        tmp<scalarField> correctionFactors() const;

public:

    //- Runtime type information
    TypeName("regionCoupled");


    // Constructors

        //- Construct from patch and internal field
        regionCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regionCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given regionCoupledFvPatchScalarField
        // onto a new patch
        regionCoupledFvPatchScalarField
        (
            const regionCoupledFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

#ifdef OPENFOAMESI
         //- Construct and return a clone
    virtual tmp<fvPatchField<scalar>> clone() const
    {
        return tmp<fvPatchField<scalar>>
        (
            new regionCoupledFvPatchScalarField(*this)
        );
    }
#endif         

        //- Construct as copy setting internal field reference
        regionCoupledFvPatchScalarField
        (
            const regionCoupledFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar> > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar> >
            (
                new regionCoupledFvPatchScalarField(*this, iF)
            );
        }     


    //- Destructor
    virtual ~regionCoupledFvPatchScalarField()
    {}


    // Member functions

        // Access
        
            //- Return true if this patch field is coupled
            virtual bool coupled() const
            {
                return coupled_;
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
                const fvPatchScalarField&,
                const labelList&
            );
            
        // Evaluation functions

            //- Return neighbour coupled internal cell data
            virtual tmp<scalarField> patchNeighbourField() const;

            //- Retrieve neighbour coupled data
            virtual void patchNeighbourField(UList<scalar>&) const;

            //- Update the coefficients associated with the patch field
            //  Sets Updated to true
            virtual void updateCoeffs();
            
            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType
            );

            //- Return patch-normal gradient
            virtual tmp<scalarField> snGrad() const;

            //- Return patch-normal gradient
            //  Note: the deltaCoeffs supplied are not used
            virtual tmp<scalarField> snGrad
            (
                const scalarField& deltaCoeffs
            ) const;
            
            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<scalarField> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<scalarField> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<scalarField> gradientInternalCoeffs
            (
                const scalarField& deltaCoeffs
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<scalarField> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<scalarField> gradientBoundaryCoeffs
            (
                const scalarField& deltaCoeffs
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<scalarField> gradientBoundaryCoeffs() const;
            
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
                Field<scalar>&,
                const Field<scalar>&,
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
                 Field<scalar>&,
                 const bool add,
                 const lduAddressing& lduAddr,
                 const label patchId,
                 const Field<scalar>&,
                 const scalarField&,
                 const Pstream::commsTypes commsType
             ) const;            
#endif            

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
