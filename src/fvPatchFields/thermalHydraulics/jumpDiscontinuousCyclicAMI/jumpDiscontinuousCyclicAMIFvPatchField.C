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

\*---------------------------------------------------------------------------*/

#include "jumpDiscontinuousCyclicAMIFvPatchField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::jumpDiscontinuousCyclicAMIFvPatchField<Type>::jumpDiscontinuousCyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    discontinuousCyclicAMIFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::jumpDiscontinuousCyclicAMIFvPatchField<Type>::jumpDiscontinuousCyclicAMIFvPatchField
(
    const jumpDiscontinuousCyclicAMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    discontinuousCyclicAMIFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::jumpDiscontinuousCyclicAMIFvPatchField<Type>::jumpDiscontinuousCyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    discontinuousCyclicAMIFvPatchField<Type>(p, iF, dict)
{
    // Call this evaluation in derived classes
    //this->evaluate(Pstream::commsTypes::buffered);
}


template<class Type>
Foam::jumpDiscontinuousCyclicAMIFvPatchField<Type>::jumpDiscontinuousCyclicAMIFvPatchField
(
    const jumpDiscontinuousCyclicAMIFvPatchField<Type>& ptf
)
:
    discontinuousCyclicAMIFvPatchField<Type>(ptf)
{}


template<class Type>
Foam::jumpDiscontinuousCyclicAMIFvPatchField<Type>::jumpDiscontinuousCyclicAMIFvPatchField
(
    const jumpDiscontinuousCyclicAMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    discontinuousCyclicAMIFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::jumpDiscontinuousCyclicAMIFvPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->primitiveField();
    const labelUList& nbrFaceCells =
        this->discontinuousCyclicAMIPatch().cyclicAMIPatch().neighbPatch().faceCells();

    Field<Type> pnf(iField, nbrFaceCells);
    tmp<Field<Type>> tpnf;

    if (this->discontinuousCyclicAMIPatch().applyLowWeightCorrection())
    {
        tpnf =
            this->discontinuousCyclicAMIPatch().interpolate
            (
                pnf,
                this->patchInternalField()()
            );
    }
    else
    {
        tpnf = this->discontinuousCyclicAMIPatch().interpolate(pnf);
    }

    if (this->doTransform())
    {
        transform(tpnf.ref(), this->forwardT(), tpnf());
    }

    tmp<Field<Type>> tjf = jump();
    if (!this->discontinuousCyclicAMIPatch().owner())
    {
        tjf = -tjf;
    }

    return tpnf - tjf;
}


template<class Type>
void Foam::jumpDiscontinuousCyclicAMIFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    // Bypass discontinuousCyclicAMI to avoid caching
    coupledFvPatchField<Type>::initEvaluate(commsType);
}


template<class Type>
void Foam::jumpDiscontinuousCyclicAMIFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Bypass discontinuousCyclicAMI to avoid caching
    coupledFvPatchField<Type>::evaluate(commsType);
    // if(this->discontinuousCyclicAMIPatch().duplicate())
    //     *this *=0.5;
}


template<class Type>
void Foam::jumpDiscontinuousCyclicAMIFvPatchField<Type>::initInterfaceMatrixUpdate
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    // Bypass discontinuousCyclicAMI to avoid caching
    coupledFvPatchField<Type>::initInterfaceMatrixUpdate
    (
        result,
        add,
        lduAddr,
        patchId,
        psiInternal,
        coeffs,
        cmpt,
        commsType
    );
}


template<class Type>
void Foam::jumpDiscontinuousCyclicAMIFvPatchField<Type>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    NotImplemented;
}


template<class Type>
void Foam::jumpDiscontinuousCyclicAMIFvPatchField<Type>::initInterfaceMatrixUpdate
(
    Field<Type>& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes commsType
) const
{
    // Bypass discontinuousCyclicAMI to avoid caching
    coupledFvPatchField<Type>::initInterfaceMatrixUpdate
    (
        result,
        add,
        lduAddr,
        patchId,
        psiInternal,
        coeffs,
        commsType
    );
}


template<class Type>
void Foam::jumpDiscontinuousCyclicAMIFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{
    const labelUList& nbrFaceCells =
        lduAddr.patchAddr
        (
            this->discontinuousCyclicAMIPatch().neighbPatchID()
        );

    Field<Type> pnf(psiInternal, nbrFaceCells);

    if (this->discontinuousCyclicAMIPatch().applyLowWeightCorrection())
    {
        pnf =
            this->discontinuousCyclicAMIPatch().interpolate
            (
                pnf,
                this->patchInternalField()()
            );

    }
    else
    {
        pnf = this->discontinuousCyclicAMIPatch().interpolate(pnf);
    }

    // only apply jump to original field
    if (&psiInternal == &this->primitiveField())
    {
        Field<Type> jf(this->jump());
        if (!this->discontinuousCyclicAMIPatch().owner())
        {
            jf *= -1.0;
        }

        pnf -= jf;
    }

    // Transform according to the transformation tensors
    this->transformCoupleField(pnf);

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add,faceCells, coeffs, pnf);
}


// ************************************************************************* //
