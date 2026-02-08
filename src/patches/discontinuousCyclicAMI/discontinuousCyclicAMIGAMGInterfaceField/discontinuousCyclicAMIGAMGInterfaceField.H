/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::discontinuousCyclicAMIGAMGInterfaceField

Description
    GAMG agglomerated cyclic interface field.

SourceFiles
    discontinuousCyclicAMIGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef discontinuousCyclicAMIGAMGInterfaceField_H
#define discontinuousCyclicAMIGAMGInterfaceField_H

#include "GAMGInterfaceField.H"
#include "discontinuousCyclicAMIGAMGInterface.H"
#include "discontinuousCyclicAMILduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class discontinuousCyclicAMIGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class discontinuousCyclicAMIGAMGInterfaceField
:
    public GAMGInterfaceField,
    virtual public discontinuousCyclicAMILduInterfaceField
{
    // Private data

        //- Local reference cast into the cyclic interface
        const discontinuousCyclicAMIGAMGInterface& discontinuousCyclicAMIInterface_;

        //- Is the transform required
        bool doTransform_;

        //- Rank of component for transformation
        int rank_;


        // Sending and receiving (distributed AMI)

            //- Current range of send requests (non-blocking)
            mutable labelRange sendRequests_;

            //- Current range of recv requests (non-blocking)
            mutable labelRange recvRequests_;

            //- Scalar send buffers
            mutable PtrList<List<solveScalar>> scalarSendBufs_;

            //- Scalar receive buffers
            mutable PtrList<List<solveScalar>> scalarRecvBufs_;


    // Private Member Functions

        //- No copy construct
        discontinuousCyclicAMIGAMGInterfaceField
        (
            const discontinuousCyclicAMIGAMGInterfaceField&
        ) = delete;

        //- No copy assignment
        void operator=(const discontinuousCyclicAMIGAMGInterfaceField&) = delete;


public:

    //- Runtime type information
    TypeName("discontinuousCyclicAMI");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        discontinuousCyclicAMIGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterfaceField
        );

        //- Construct from GAMG interface and fine level interface field
        discontinuousCyclicAMIGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const bool doTransform,
            const int rank
        );

        //- Construct from GAMG interface and Istream
        discontinuousCyclicAMIGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            Istream& is
        );

        //- Construct from GAMG interface and local and remote fields
        discontinuousCyclicAMIGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& local,
            const UPtrList<lduInterfaceField>& other
        );

        //- Construct by assembling and return a clone.
        virtual autoPtr<GAMGInterfaceField> clone
        (
            const GAMGInterface& GAMGCp,
            const UPtrList<lduInterfaceField>& other
        ) const
        {
            return autoPtr<GAMGInterfaceField>
            (
                new discontinuousCyclicAMIGAMGInterfaceField
                (
                    GAMGCp,
                    *this,      // local field
                    other       // other fields
                )
            );
        }


    //- Destructor
    virtual ~discontinuousCyclicAMIGAMGInterfaceField() = default;


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return discontinuousCyclicAMIInterface_.size();
            }


        // Interface matrix update

            //- Are all (receive) data available?
            virtual bool ready() const;

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


        //- Cyclic interface functions

            //- Does the interface field perform the transformation
            virtual bool doTransform() const
            {
                return doTransform_;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return discontinuousCyclicAMIInterface_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return discontinuousCyclicAMIInterface_.reverseT();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return rank_;
            }


        // I/O

            //- Write to stream
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
