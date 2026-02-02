/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::discontinuousCyclicAMIFvPatch

Description
    Cyclic patch for Arbitrary Mesh Interface (AMI)

SourceFiles
    discontinuousCyclicAMIFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_discontinuousCyclicAMIFvPatch_H
#define Foam_discontinuousCyclicAMIFvPatch_H

#include "discontinuousCyclicAMIPolyPatch.H"
#include "coupledFvPatch.H"
#include "discontinuousCyclicAMILduInterface.H"
#include "fvBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class discontinuousCyclicAMIFvPatch Declaration
\*---------------------------------------------------------------------------*/

class discontinuousCyclicAMIFvPatch
:
    public coupledFvPatch,
    public discontinuousCyclicAMILduInterface
{
    // Private Data

        const discontinuousCyclicAMIPolyPatch& discontinuousCyclicAMIPolyPatch_;

protected:

    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;

        //- Correct patch deltaCoeffs
        virtual void makeDeltaCoeffs(scalarField&) const;

        //- Correct patch non-ortho deltaCoeffs
        virtual void makeNonOrthoDeltaCoeffs(scalarField&) const;

        //- Correct patch non-ortho correction vectors
        virtual void makeNonOrthoCorrVectors(vectorField&) const;

        //- Correct patches after moving points
        virtual void movePoints();


public:

    //- Runtime type information
    TypeName(discontinuousCyclicAMIPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        discontinuousCyclicAMIFvPatch
        (
            const polyPatch& patch, 
            const fvBoundaryMesh& bm
        );



            //- Return local reference cast into the cyclic patch
            const discontinuousCyclicAMIPolyPatch& cyclicAMIPatch() const
            {
                return discontinuousCyclicAMIPolyPatch_;
            }

            //- Does this side own the patch?
            virtual bool owner() const
            {
                return discontinuousCyclicAMIPolyPatch_.owner();
            }

            //- Return neighbour patch ID
            virtual label neighbPatchID() const
            {
                return discontinuousCyclicAMIPolyPatch_.neighbPatchID();
            }

            //- Return a reference to the neighbour patch
            virtual const discontinuousCyclicAMIFvPatch& neighbPatch() const
            {
                return refCast<const discontinuousCyclicAMIFvPatch>
                (
                    this->boundaryMesh()[discontinuousCyclicAMIPolyPatch_.neighbPatchID()]
                );
            }

            //- Return a reference to the AMI interpolator
            virtual const AMIPatchToPatchInterpolation& AMI() const
            {
                return discontinuousCyclicAMIPolyPatch_.AMI();
            }

            //- Return true if applying the low weight correction
            virtual bool applyLowWeightCorrection() const
            {
                return discontinuousCyclicAMIPolyPatch_.applyLowWeightCorrection();
            }


            //- Are the cyclic planes parallel
            virtual bool parallel() const
            {
                return discontinuousCyclicAMIPolyPatch_.parallel();
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return discontinuousCyclicAMIPolyPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return discontinuousCyclicAMIPolyPatch_.reverseT();
            }

            const discontinuousCyclicAMIFvPatch& neighbFvPatch() const
            {
                return refCast<const discontinuousCyclicAMIFvPatch>
                (
                    this->boundaryMesh()[discontinuousCyclicAMIPolyPatch_.neighbPatchID()]
                );
            }


            //- Return true if this patch is coupled. This is equivalent
            //- to the coupledPolyPatch::coupled() if parallel running or
            //- both sides present, false otherwise
            virtual bool coupled() const;

            virtual bool duplicate() const
            {
                return discontinuousCyclicAMIPolyPatch_.duplicate();
            }

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;

            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const Field<Type>& fld,
                const UList<Type>& defaultValues = UList<Type>()
            ) const
            {
                return discontinuousCyclicAMIPolyPatch_.interpolate(fld, defaultValues);
            }

            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const tmp<Field<Type>>& tFld,
                const UList<Type>& defaultValues = UList<Type>()
            ) const
            {
                return discontinuousCyclicAMIPolyPatch_.interpolate(tFld, defaultValues);
            }


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //- the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;


            //- Return the values of the given internal data adjacent to
            //- the interface as a field using a mapping faceCell
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData,
                const labelUList& faceCells
            ) const;


            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
