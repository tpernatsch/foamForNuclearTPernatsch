/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "timeDependentAverageMappedFvPatchField.H"
#include "volFields.H"
#include "interpolationCell.H"
#include "pTraits.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::timeDependentAverageMappedFvPatchField<Type>
    ::timeDependentAverageMappedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    mappedPatchFieldBase<Type>
    (
        this->mapper(p, iF),
        *this
    ),
    timePoints_(),
    averageValues_(),
    method_(interpolationMethod::LINEAR),
    table_(timePoints_, averageValues_, method_)
{}


template<class Type>
Foam::timeDependentAverageMappedFvPatchField<Type>
    ::timeDependentAverageMappedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    mappedPatchFieldBase<Type>
    (
        this->mapper(p, iF),
        *this,
        dict
    ),
    timePoints_(dict.lookup("timePoints")),
    averageValues_(dict.lookup("averageValues")),
    method_
    (
        interpolateTableBase::interpolationMethodNames_
        [
            dict.lookupOrDefault<word>
            ("timeInterpolationMethod", "linear")
        ]
    ),
    table_(timePoints_, averageValues_, method_)
{}


template<class Type>
Foam::timeDependentAverageMappedFvPatchField<Type>
    ::timeDependentAverageMappedFvPatchField
(
    const timeDependentAverageMappedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    mappedPatchFieldBase<Type>(this->mapper(p, iF), *this, ptf),
    timePoints_(ptf.timePoints_),
    averageValues_(ptf.averageValues_),
    method_(ptf.method_),
    table_(timePoints_, averageValues_, method_)
{}


template<class Type>
Foam::timeDependentAverageMappedFvPatchField<Type>
    ::timeDependentAverageMappedFvPatchField
(
    const timeDependentAverageMappedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    mappedPatchFieldBase<Type>(this->mapper(this->patch(), iF), *this, ptf),
    timePoints_(ptf.timePoints_),
    averageValues_(ptf.averageValues_),
    method_(ptf.method_),
    table_(timePoints_, averageValues_, method_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::timeDependentAverageMappedFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
}


template<class Type>
void Foam::timeDependentAverageMappedFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);
}


template<class Type>
const Foam::mappedPatchBase& Foam
    ::timeDependentAverageMappedFvPatchField<Type>::mapper
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
{
    if (!isA<mappedPatchBase>(p.patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.patch().name()
            << " of field " << iF.name()
            << " in file " << iF.objectPath()
            << exit(FatalError);
    }
    return refCast<const mappedPatchBase>(p.patch());
}


template<class Type>
void Foam::timeDependentAverageMappedFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar time = this->db().time().value();
    const scalar userTime = this->db().time().timeToUserTime(time);
    Type average = this->table_(userTime);
    
    // Because mappedPatchFieldBase declares these as constant, we must
    // force the changes to these variables    
    *((bool*)const_cast<bool*>(&this->setAverage_)) = true;
    *((Type*)const_cast<Type*>(&this->average_)) = average;

    this->operator==(this->mappedField());

    if (debug)
    {
        Info<< "mapped on field:"
            << this->internalField().name()
            << " patch:" << this->patch().name()
            << "  avg:" << gAverage(*this)
            << "  min:" << gMin(*this)
            << "  max:" << gMax(*this)
            << endl;
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void
Foam::timeDependentAverageMappedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    mappedPatchFieldBase<Type>::write(os);

#ifdef OPENFOAMFOUNDATION
    writeEntry(os, "value", *this);
    writeEntry(os, "timePoints", 
               static_cast<const List<scalar>&>(timePoints_));
    writeEntry(os, "averageValues", 
               static_cast<const List<Type>&>(averageValues_));
    writeEntry(os, "timeInterpolationMethod", 
               interpolateTableBase::interpolationMethodNames_[method_]);    
#elif OPENFOAMESI
    this->writeEntry("value", os);
    os.writeEntry<List<scalar>>("timePoints", 
            static_cast<const List<scalar>&>(timePoints_));
    os.writeEntry<List<Type>>("averageValues", 
            static_cast<const List<Type>&>(averageValues_));
    os.writeEntry<word>("timeInterpolationMethod", 
            interpolateTableBase::interpolationMethodNames_[method_]);
#endif
}


// ************************************************************************* //
