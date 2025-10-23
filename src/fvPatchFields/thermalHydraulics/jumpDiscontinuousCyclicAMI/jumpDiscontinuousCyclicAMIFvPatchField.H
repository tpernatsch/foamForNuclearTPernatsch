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
    Foam::jumpDiscontinuousCyclicAMIFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition provides a base class that enforces a cyclic
    condition with a specified 'jump' (or offset) between a pair of boundaries,
    whereby communication between the patches is performed using an arbitrary
    mesh interface (AMI) interpolation. Identical to jumpCyclicAMIFvPatchField,
    but for discontinuous cyclic AMI patches.

See also
    Foam::discontinuousCyclicAMIFvPatchField

SourceFiles
    jumpDiscontinuousCyclicAMIFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef jumpDiscontinuousCyclicAMIFvPatchField_H
#define jumpDiscontinuousCyclicAMIFvPatchField_H

#include "discontinuousCyclicAMIFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class jumpDiscontinuousCyclicAMIFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class jumpDiscontinuousCyclicAMIFvPatchField
:
    public discontinuousCyclicAMIFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("jumpdiscontinuousCyclicAMI");


    // Constructors

        //- Construct from patch and internal field
        jumpDiscontinuousCyclicAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        jumpDiscontinuousCyclicAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given jumpDiscontinuousCyclicAMIFvPatchField onto a
        //  new patch
        jumpDiscontinuousCyclicAMIFvPatchField
        (
            const jumpDiscontinuousCyclicAMIFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        jumpDiscontinuousCyclicAMIFvPatchField
        (
            const jumpDiscontinuousCyclicAMIFvPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        jumpDiscontinuousCyclicAMIFvPatchField
        (
            const jumpDiscontinuousCyclicAMIFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );


    // Member functions

        // Access

            virtual bool coupled() const
            {
                return true;
            }

            //- Return the interface type
            virtual const word& interfaceFieldType() const
            {
                return discontinuousCyclicAMIFvPatchField<Type>::type();
            }

            //- Return the "jump" across the patch as a "half" field
            virtual tmp<Field<Type>> jump() const = 0;


        // Evaluation functions

            //- Return neighbour coupled given internal cell data
            virtual tmp<Field<Type>> patchNeighbourField() const;

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

            //- Evaluate the patch field
            virtual void evaluate(const Pstream::commsTypes commsType);

            //- Initialise interface functionality
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

            //- Initialise interface functionality
            void initInterfaceMatrixUpdate
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
};


//- Update result field based on interface functionality
template<>
void jumpDiscontinuousCyclicAMIFvPatchField<scalar>::updateInterfaceMatrix
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "jumpDiscontinuousCyclicAMIFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
