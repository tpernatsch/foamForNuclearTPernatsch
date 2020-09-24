/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "patchScalarFieldValue.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(patchScalarFieldValue, 0);
    addToRunTimeSelectionTable
    (
        functionObject, 
        patchScalarFieldValue, 
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::patchScalarFieldValue::writeFileHeader(Ostream& os)
{
    /*
    if (writtenHeader_)
    {
        os << endl;
    }
    else
    {
        word headerText
        (
            "Value of volScalarField "+fieldName_+" on patch "+patchName_
        );
        writeHeader(os, headerText);
    }
    
    word time("Time = "+mesh_.time().timeName());

    os << time << endl;
    */

    writtenHeader_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::patchScalarFieldValue::patchScalarFieldValue
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    fieldPtr_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::patchScalarFieldValue::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict) && writeFile::read(dict))
    {
        fieldName_ = dict.get<word>("field");

        patchName_ = dict.get<word>("patch");

        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        
        patchID_ = pbm.findIndex(patchName_);
        
        return true;
    }

    return false;
}


bool Foam::functionObjects::patchScalarFieldValue::execute()
{
    return true;
}


bool Foam::functionObjects::patchScalarFieldValue::write()
{
    //writeFileHeader(file());

    file() << endl;

    if (fieldPtr_ == nullptr)
    {
        fieldPtr_ = &(mesh_.lookupObject<volScalarField>(fieldName_));
    }
    const volScalarField& field(*fieldPtr_);
    const scalarField& bField(field.boundaryField()[patchID_]);

    const scalar t = mesh_.time().timeOutputValue();

    file() << "( " << t << " ( "; 
    forAll(bField, i)
    {
        file() << bField[i] << " ";
    }
    file() << "))";
    
    return true;
}


// ************************************************************************* //
