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

/*---------------------------------------------------------------------------*\
Class
    Foam::connectingJunctionJumpAMIFvPatchVectorField

Group
    grpCoupledBoundaryConditions

Author
    Giovanni Nervi, EPFL (giovanni.nervi@epfl.ch)

Description
    This boundary condition handles the connection of multiple 1D branches
    (e.g. pipes) into a single outlet across non-conformal AMI interfaces.
    It is the counterpart of branchingJunctionJumpAMI, merging an arbitrary
    number of incoming flows into one patch.

    Rather than prescribing a fixed jump, the condition balances flow and
    pressure among the connected patches and ensures consistent coupling
    at the junction.

Usage
    \table
        Property         | Description                                                    | Required | Default
        patchType        | underlying patch type (should be 'discontinuousCyclicAMI')     | yes      |
        connectedPatches | list of patches representing the other branches connected here | yes      |
        underRelaxation  | relaxation factor for iterative coupling                       | no       | 1
    \endtable

    Example:
    \verbatim
    <patchName>
    {
        type              connectingJunctionJumpAMI;
        patchType         discontinuousCyclicAMI;
        connectedPatches  (pipe2AMI pipe3AMI pipe4AMI);
        underRelaxation   0.5;
    }
    \endverbatim

Note
    The underlying patchType must be 'discontinuousCyclicAMI'. The
    connectedPatches entry should list all other branches that join into
    this patch, excluding the neighbour already paired by the AMI.

See also
    Foam::jumpDiscontinuousCyclicAMIFvPatchField

SourceFiles
    connectingJunctionJumpAMIFvPatchVectorField.C
\*---------------------------------------------------------------------------*/


#ifndef connectingJunctionJumpAMIFvPatchVectorField_H
#define connectingJunctionJumpAMIFvPatchVectorField_H

#include "jumpDiscontinuousCyclicAMIFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class connectingJunctionJumpAMIFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/


class connectingJunctionJumpAMIFvPatchVectorField
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
    TypeName("connectingJunctionJumpAMI");


    // Constructors

        //- Construct from patch and internal field
        connectingJunctionJumpAMIFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        connectingJunctionJumpAMIFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given connectingJunctionJumpAMIFvPatchVectorField onto a
        //  new patch
        connectingJunctionJumpAMIFvPatchVectorField
        (
            const connectingJunctionJumpAMIFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        connectingJunctionJumpAMIFvPatchVectorField
        (
            const connectingJunctionJumpAMIFvPatchVectorField&
        );

        //- Construct as copy setting internal field reference
        connectingJunctionJumpAMIFvPatchVectorField
        (
            const connectingJunctionJumpAMIFvPatchVectorField&,
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
//     #include "connectingJunctionJumpAMIFvPatchVectorField.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
