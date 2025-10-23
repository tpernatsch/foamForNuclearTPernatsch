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
    Foam::branchingJunctionJumpAMIFvPatchVectorField

Group
    grpCoupledBoundaryConditions

Author
    Giovanni Nervi, EPFL (giovanni.nervi@epfl.ch)

Description
    This boundary condition handles flow branching across non-conformal
    cyclic AMI interfaces. It is intended for 1D (pipe-type) networks where
    one incoming branch splits into an arbitrary number of outgoing branches.

    Instead of prescribing a fixed jump, this condition distributes and
    balances the flow and pressure across the connected patches according
    to the specified topology.

    The user specifies:
      - the underlying patch type (as for other discontinuous cyclic patches),
      - the list of other branch patches connected to this one,
      - optional under-relaxation for the iterative update.

Usage
    \table
        Property         | Description                                                    | Required | Default
        patchType        | underlying patch type (should be 'discontinuousCyclicAMI')     | yes      |
        conectedPatches  | list of patches representing the other branches                | yes      |
        underRelaxation  | relaxation factor for iterative coupling                       | no       | 1
    \endtable

    Example:
    \verbatim
    <patchName>
    {
        type            branchingJunctionJumpAMI;
        patchType       discontinuousCyclicAMI;
        conectedPatches   (pipe2AMI pipe3AMI pipe4AMI); // all other patches except for nbrPatch
        underRelaxation 0.5;
    }
    \endverbatim

Note
    The underlying patchType must be 'discontinuousCyclicAMI'. The
    conectedPatches entry should list all the additional branches that form
    part of the junction, excluding the neighbour patch directly paired by
    the AMI.

See also
    Foam::jumpDiscontinuousCyclicAMIFvPatchField

SourceFiles
    branchingJunctionJumpAMIFvPatchVectorField.C
\*---------------------------------------------------------------------------*/

#ifndef branchingJunctionJumpAMIFvPatchVectorField_H
#define branchingJunctionJumpAMIFvPatchVectorField_H

#include "jumpDiscontinuousCyclicAMIFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class branchingJunctionJumpAMIFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/


class branchingJunctionJumpAMIFvPatchVectorField
:
    public jumpDiscontinuousCyclicAMIFvPatchField<vector>
{

protected:

    // Protected data

    vectorField jump_;

    scalar underRelaxation_;

    scalarList scalingFactors_;

    wordList overlapPatchNames_;

public:

    //- Runtime type information
    TypeName("branchingJunctionJumpAMI");


    // Constructors

        //- Construct from patch and internal field
        branchingJunctionJumpAMIFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        branchingJunctionJumpAMIFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given branchingJunctionJumpAMIFvPatchVectorField onto a
        //  new patch
        branchingJunctionJumpAMIFvPatchVectorField
        (
            const branchingJunctionJumpAMIFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        branchingJunctionJumpAMIFvPatchVectorField
        (
            const branchingJunctionJumpAMIFvPatchVectorField&
        );

        //- Construct as copy setting internal field reference
        branchingJunctionJumpAMIFvPatchVectorField
        (
            const branchingJunctionJumpAMIFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return fvPatchField<vector>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return fvPatchField<vector>::Clone(*this, iF);
        }

        //         //- Return neighbour coupled given internal cell data
        // virtual tmp<Field<vector>> patchNeighbourField() const;

        //- Update result field based on interface functionality
        // virtual void updateInterfaceMatrix
        // (
        //     Field<vector>&,
        //     const bool add,
        //     const lduAddressing& lduAddr,
        //     const label patchId,
        //     const Field<vector>&,
        //     const scalarField&,
        //     const Pstream::commsTypes commsType
        // ) const;



    // Member functions

        // Access

            //- Return the "jump" across the patch
            virtual tmp<Field<vector>> jump() const;


        // Mapping functions
            // virtual void evaluate
            // (
            //     const Pstream::commsTypes commsType
            // );
            virtual void updateCoeffs();

            scalarList scalingFactors() const
            {
                return scalingFactors_;
            }

            wordList overlappingPatch() const
            {
                return overlapPatchNames_;
            }

            // virtual void evaluate(const Pstream::commsTypes);


            virtual void updateInterfaceMatrix
            (
                solveScalarField& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes
            ) const;


            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap(const fvPatchField<vector>&, const labelList&);


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
//     #include "branchingJunctionJumpAMIFvPatchVectorField.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
