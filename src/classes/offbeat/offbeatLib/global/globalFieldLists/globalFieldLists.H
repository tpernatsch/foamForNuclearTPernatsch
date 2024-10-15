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

Global
    globalFieldLists

Description
    PtrLists of scalar/vector/symmTensor fields that can be accessed from 
    registry.

    The relax, storePrevIter and correctBoundaryConditions can be called 
    on all the registered fields.

SourceFiles
    globalFieldLists.C

\*---------------------------------------------------------------------------*/

#ifndef globalFieldLists_H
#define globalFieldLists_H

#include "volFields.H"
#include "dimensionSet.H"
#include "fvMesh.H"
#include "calculatedFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
//- List of volumetric fields
PtrList<Foam::GeometricField<Type, Foam::fvPatchField, 
        Foam::volMesh>> volFields_;

//- Create field or lookup from registry and return non-const reference
template<class Type>
Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& 
createOrLookup
(
    const fvMesh& mesh,
    const word& name,
    const dimensionSet& dimensions = dimless,
    const Type initialValue = pTraits<Type>::zero,
    const word & patchFieldType = fvPatchField<Type>::calculatedType()
)
{    
    if(mesh.foundObject<GeometricField<Type, Foam::fvPatchField, 
        Foam::volMesh>>(name))
    {
        return 
        mesh.lookupObjectRef<Foam::GeometricField<Type, Foam::fvPatchField, 
        Foam::volMesh>>(name);
    }
    else
    {
        volFields_<Type>.append
        (
            new Foam::GeometricField<Type, fvPatchField, volMesh>
            (
                IOobject
                (
                    name,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensioned<Type>(name, dimensions, initialValue),
                patchFieldType
            )
        );

        return volFields_<Type>[volFields_<Type>.size() - 1];

    }
};

template<class Type>
//- List of fields
PtrList<Foam::Field<Type>> fields_;

//- Create field or read from dictionary
template<class Type>
Foam::Field<Type>& 
createOrRead
(
    const fvMesh& mesh,
    const word& name,
    const dictionary& dict,
    const bool& addToRegistry, 
    const Type initialValue = pTraits<Type>::zero
)
{    
    //- If addToRegistry is true, create geometric field and read internal
    //  field
    if(addToRegistry)
    {
        Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>&
        geometricField
        (
            createOrLookup<Type>
            (
                mesh, name, dimless, initialValue, 
                zeroGradientFvPatchField<Type>::typeName
            )
        );

        return geometricField.ref();
    }

    // Read the field from geometric field file if present
    IOobject geometricFieldIO
    (
        name,
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE,
        false
    );

    //- If the field is present in the folder, read the internal field
    if
    (
        geometricFieldIO.template 
        typeHeaderOk<GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
        (
            true
        )
    )
    {
        autoPtr<GeometricField<Type, Foam::fvPatchField, Foam::volMesh>> 
        geometricField
        (
            new GeometricField<Type, Foam::fvPatchField, Foam::volMesh>
            (
                geometricFieldIO,
                mesh,
                dimensioned<Type>(name, dimless, initialValue)
            )
        );

        fields_<Type>.append
        (
            new Foam::Field<Type>(geometricField->internalField())
        );

        return fields_<Type>[fields_<Type>.size() - 1];
    }
    
    //- Read, if the field is present in the dictionary
#ifdef OPENFOAMFOUNDATION
    if(!dict.isNull() and dict.found(name))
#elif OPENFOAMESI
    if(!dict.isNullDict() and dict.found(name))
#endif                
    {
        fields_<Type>.append
        (
            new Foam::Field<Type>(dict.lookup(name))
        );

        if(fields_<Type>[fields_<Type>.size() - 1].size() != mesh.nCells())
        {             
            FatalErrorIn("fgrSCIANTIX::readField(const word&, const Type")
            << "Size of field " << name << " read from regIOobject "
            << "is different than the total number of cells." 
            << exit(FatalError);
        }

        return fields_<Type>[fields_<Type>.size() - 1];
    }
    //- Create, if the field is not present in the dictionary
#ifdef OPENFOAMESIORFOUNDATION
    else
#endif 
    {
        fields_<Type>.append
        (
            new Foam::Field<Type>(mesh.nCells(), initialValue)
        );

        return fields_<Type>[fields_<Type>.size() - 1];
    }
};

//- Store previous iteration of global fields
//- TODO: check if it is necessary as it forces the code to create prevIter
//- fields
template<class Type>
void storePrevIter()
{
    forAll(volFields_<Type>, i)
    {
        volFields_<Type>[i].storePrevIter();
    }        
};

//- Call storePrevIter() for each GeometricField<Type> PtrList
void storeGlobalFieldsPrevIter();

//- Relax global fields
//- TODO: check if it is necessary as it is done for every field
template<class Type>
void relax()
{
    forAll(volFields_<Type>, i)
    {
        volFields_<Type>[i].relax();
    }        
};

//- Call relax() for each GeometricField<Type> PtrList
void relaxGlobalFields();

//- Correct boundary conditions
//- TODO: check if it is necessary as it is done for every field
template<class Type>
void correctBoundaryConditions()
{
    forAll(volFields_<Type>, i)
    {
        volFields_<Type>[i].correctBoundaryConditions();
    }        
};

//- Call relax() for each GeometricField<Type> PtrList
void correctBCsGlobalFields();


} // End namespace Foam
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
//  #include "globalFieldLists.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif