/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#ifndef sliceMapperTemplates_H
#define sliceMapperTemplates_H

#include "sliceMapper.H"
#include "fvc.H"
#include "reverseLinear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type>
//- List of average fields for template type
autoPtr<Foam::HashTable<Foam::Field<Type>, Foam::word>> avgFieldsPtr_;

template<class Type>
const Foam::Field<Type>&
Foam::sliceMapper::sliceAverage
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{   
    // If not set yet, create HashTable of average fields for template type
    if(!avgFieldsPtr_<Type>.valid())
    {
        avgFieldsPtr_<Type>.reset
        (
            new 
            HashTable<Foam::Field<Type>, word>(10)
        );
    }

    // Reference to HashTable of average fields for template type
    HashTable<Foam::Field<Type>, word>& avgFields = avgFieldsPtr_<Type>();

    // Name of volumetric field passed in argument
    word fieldName(vf.name());

    // If vf name is not present in HashTable, insert a field with all zeroes
    if(avgFields.find(fieldName) == avgFields.end())
    {
        avgFields.insert
        (
            fieldName, 
            Foam::Field<Type>(nSlices(), pTraits<Type>::zero)
        );
    }

    // Internal field reference of vf passed by argument
    const Field<Type>& vfI = vf.internalField();

    // Reference to average field corresponding to vf name
    Foam::Field<Type>& avgField(avgFields[fieldName]);

    // Reference to mesh volumes
    const scalarField& Vi = mesh_.V();

    const PtrList<labelList>& sliceAddr(sliceAddrList());

    // Loop over slices
    forAll(sliceAddr, i)
    {
        // Addr list for current slice
        const labelList& sliceAddri = sliceAddr[i];

        // Map vf internal field in current slice
        Foam::Field<Type> sliceField(vfI, sliceAddri);
        scalarField sliceV(Vi, sliceAddri);

        // Set volumetric average for current slice
        avgField[i] = gSum(sliceField*sliceV)/gSum(sliceV);
    }

    return avgField;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

