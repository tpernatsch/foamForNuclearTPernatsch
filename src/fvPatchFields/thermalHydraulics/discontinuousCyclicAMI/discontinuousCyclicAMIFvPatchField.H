/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2506                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2025 OpenCFD Ltd.          |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::discontinuousCyclicAMIFvPatchField

Group
    grpCoupledBoundaryConditions

Author
    Giovanni Nervi (giovanni.nervi@epfl.ch)

Description
    This boundary condition is a modification of the OpenFOAM cyclicAMIFvPatchField.
    It switches between a zeroGradient BC and a jump cyclicBC, depending on whether 
    the patchField is being created from dictionary (i.e., read from file) or being 
    created during the simulation. The two modes of the patchField allow to have 
    discontinuity in certain fields (e.g., porosity, hydraulic diameter)
    which are not physically continuous across AMI interface, while still allowing
    other fields (e.g., velocity, pressure) to be coupled.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            discontinuousCyclicAMI;
    }
    \endverbatim

Note
    The outer boundary of the patch pairs must be similar, i.e. if the owner
    patch is transformed to the neighbour patch, the outer perimeter of each
    patch should be identical (or very similar).

    The \c neighbourValue is only needed when running distributed,
    i.e. the neighbour cells are on a different processor from the owner cells.
    It guarantees consistent restart e.g. when doing a snGrad and avoids
    additional communication.

See also
    Foam::AMIInterpolation
    Foam::processorFvPatchField

SourceFiles
    discontinuousCyclicAMIFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_discontinuousCyclicAMIFvPatchField_H
#define Foam_discontinuousCyclicAMIFvPatchField_H

#include "discontinuousCyclicAMIFvPatch.H"
#include "coupledFvPatchField.H"
#include "discontinuousCyclicAMILduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class discontinuousCyclicAMIFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class discontinuousCyclicAMIFvPatchField
:
    virtual public discontinuousCyclicAMILduInterfaceField,
    public coupledFvPatchField<Type>
{
    // Private Data

        //- Local reference cast into the cyclic patch
        const discontinuousCyclicAMIFvPatch& discontinuousCyclicAMIPatch_;


        // Sending and receiving (distributed AMI)

            //- Current range of send requests (non-blocking)
            mutable labelRange sendRequests_;

            //- Current range of recv requests (non-blocking)
            mutable labelRange recvRequests_;

            //- Send buffers
            mutable PtrList<List<Type>> sendBufs_;

            //- Receive buffers
            mutable PtrList<List<Type>> recvBufs_;

            //- Scalar send buffers
            mutable PtrList<List<solveScalar>> scalarSendBufs_;

            //- Scalar receive buffers
            mutable PtrList<List<solveScalar>> scalarRecvBufs_;


            // Only used for AMI caching

                //- Current range of send requests (non-blocking)
                mutable labelRange sendRequests1_;

                //- Current range of recv requests (non-blocking)
                mutable labelRange recvRequests1_;

                //- Send buffers
                mutable PtrList<List<Type>> sendBufs1_;

                //- Receive buffers_
                mutable PtrList<List<Type>> recvBufs1_;

                //- Scalar send buffers
                mutable PtrList<List<solveScalar>> scalarSendBufs1_;

                //- Scalar receive buffers
                mutable PtrList<List<solveScalar>> scalarRecvBufs1_;

            //- Neighbour coupled internal cell data
            mutable autoPtr<Field<Type>> patchNeighbourFieldPtr_;

            bool coupled_;

            bool applyDuplicate_;



    // Private Member Functions

        //- Return the AMI corresponding to the owner side
        const AMIPatchToPatchInterpolation& ownerAMI() const
        {
            return
            (
                discontinuousCyclicAMIPatch_.owner()
              ? discontinuousCyclicAMIPatch_.AMI()
              : discontinuousCyclicAMIPatch_.neighbPatch().AMI()
            );
        }

        //- All receive/send requests have completed
        virtual bool all_ready() const;

        //- Use neighbour field caching
        static bool cacheNeighbourField();

        //- Return neighbour coupled internal cell data
        virtual tmp<Field<Type>> getNeighbourField(const UList<Type>&) const;

        //- Return neighbour coupled internal cell data (cached or extracted),
        //- with optional check that AMI.comm() is valid.
        tmp<Field<Type>> getPatchNeighbourField(bool checkCommunicator) const;

        //- Return new matrix coeffs
        tmp<Field<scalar>> coeffs
        (
            fvMatrix<Type>& matrix,
            const Field<scalar>&,
            const label
        ) const;

        template<class Type2>
        void collectStencilData
        (
            const refPtr<mapDistribute>& mapPtr,
            const labelListList& stencil,
            const Type2& data,
            List<Type2>& expandedData
        );


public:

    //- Runtime type information
    TypeName(discontinuousCyclicAMIFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        discontinuousCyclicAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        discontinuousCyclicAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given discontinuousCyclicAMIFvPatchField onto a new patch
        discontinuousCyclicAMIFvPatchField
        (
            const discontinuousCyclicAMIFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        discontinuousCyclicAMIFvPatchField(const discontinuousCyclicAMIFvPatchField<Type>&);

        //- Construct as copy setting internal field reference
        discontinuousCyclicAMIFvPatchField
        (
            const discontinuousCyclicAMIFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return fvPatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return fvPatchField<Type>::Clone(*this, iF);
        }


    // Member Functions

        //- Return local reference cast into the cyclic AMI patch
        const discontinuousCyclicAMIFvPatch& discontinuousCyclicAMIPatch() const noexcept
        {
            return discontinuousCyclicAMIPatch_;
        }


    // Coupling

        //- Major difference w.r.t. cyclicAMI --> coupled is defaulted to false
        // unless found in dictionary. If false, then the patch defaults to a 
        // zeroGradient BC. This is done to allow to have discontinuity between
        // the default fields on the two patches connected by the discontinuous 
        // cyclic AMI

        virtual bool coupled() const
        {
            return coupled_;
        }

        virtual void coupleField()
        {
            coupled_ = true;
        }

        virtual void deCoupleField()
        {
            coupled_ = false;
        }

        //- Are all (receive) data available?
        virtual bool ready() const;

        //- Return neighbour coupled internal cell data
        virtual tmp<Field<Type>> patchNeighbourField() const;

        //- Retrieve neighbour coupled internal cell data
        virtual void patchNeighbourField(UList<Type>& pnf) const;

        //- Return reference to neighbour patchField
        const discontinuousCyclicAMIFvPatchField<Type>& neighbourPatchField() const;


    // Mapping Functions

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap
        (
            const fvPatchFieldMapper&
        );

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap
        (
            const fvPatchField<Type>&,
            const labelList&
        );


    // Evaluation

        //- Initialise the evaluation of the patch field
        virtual void initEvaluate(const Pstream::commsTypes commsType);

        //- Evaluate the patch field
        virtual void evaluate(const Pstream::commsTypes commsType);

        //- Return gradient at boundary
        virtual tmp<Field<Type>> snGrad
        (
            const scalarField& deltaCoeffs
        ) const;

        virtual tmp<Field<Type>> snGrad() const
        {
            return tmp<Field<Type>>::New(this->size(), Zero);
        }

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<Field<Type>> valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<Field<Type>> valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const;

        
        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<Type>> gradientInternalCoeffs
        (
            const scalarField& deltaCoeffs
        ) const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<Type>> gradientBoundaryCoeffs
        (
            const scalarField& deltaCoeffs
        ) const;

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<Type>> gradientInternalCoeffs() const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<Type>> gradientBoundaryCoeffs() const;

        // Coupled interface functionality

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                solveScalarField& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                solveScalarField& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                Field<Type>& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const Field<Type>& psiInternal,
                const scalarField& coeffs,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>&,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const Field<Type>&,
                const scalarField&,
                const Pstream::commsTypes commsType
            ) const;


            //- Manipulate matrix
            virtual void manipulateMatrix
            (
                fvMatrix<Type>& m,
                const label iMatrix,
                const direction cmpt
            );

        // Coupled interface functions

            //- Does the patch field perform the transformation
            virtual bool doTransform() const
            {
                return (pTraits<Type>::rank && !discontinuousCyclicAMIPatch_.parallel());
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return discontinuousCyclicAMIPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return discontinuousCyclicAMIPatch_.reverseT();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return pTraits<Type>::rank;
            }


        // I-O

            //- Write
            virtual void write(Ostream& os) const;


        // Member Operators

            virtual void operator=(const fvPatchField<Type>&);
            virtual void operator==(const fvPatchField<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "discontinuousCyclicAMIFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
