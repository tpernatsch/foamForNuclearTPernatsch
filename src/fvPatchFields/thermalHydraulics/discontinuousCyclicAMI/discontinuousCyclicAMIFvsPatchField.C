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

#include "discontinuousCyclicAMIFvsPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::discontinuousCyclicAMIFvsPatchField<Type>::discontinuousCyclicAMIFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(p, iF),
    discontinuousCyclicAMIPatch_(refCast<const discontinuousCyclicAMIFvPatch>(p))
{}


template<class Type>
Foam::discontinuousCyclicAMIFvsPatchField<Type>::discontinuousCyclicAMIFvsPatchField
(
    const discontinuousCyclicAMIFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvsPatchField<Type>(ptf, p, iF, mapper),
    discontinuousCyclicAMIPatch_(refCast<const discontinuousCyclicAMIFvPatch>(p))
{
    if (!isA<discontinuousCyclicAMIFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
Foam::discontinuousCyclicAMIFvsPatchField<Type>::discontinuousCyclicAMIFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    coupledFvsPatchField<Type>(p, iF, dict),
    discontinuousCyclicAMIPatch_(refCast<const discontinuousCyclicAMIFvPatch>(p, dict))
{
    if (!isA<discontinuousCyclicAMIFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "patch " << this->patch().index() << " not discontinuousCyclicAMI type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::discontinuousCyclicAMIFvsPatchField<Type>::discontinuousCyclicAMIFvsPatchField
(
    const discontinuousCyclicAMIFvsPatchField<Type>& ptf
)
:
    coupledFvsPatchField<Type>(ptf),
    discontinuousCyclicAMIPatch_(ptf.discontinuousCyclicAMIPatch_)
{}


template<class Type>
Foam::discontinuousCyclicAMIFvsPatchField<Type>::discontinuousCyclicAMIFvsPatchField
(
    const discontinuousCyclicAMIFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(ptf, iF),
    discontinuousCyclicAMIPatch_(ptf.discontinuousCyclicAMIPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::discontinuousCyclicAMIFvsPatchField<Type>::coupled() const
{
    return discontinuousCyclicAMIPatch_.coupled();
}


// ************************************************************************* //
