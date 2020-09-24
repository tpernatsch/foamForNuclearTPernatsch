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

#include "stopIfMaxFieldDiff.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(stopIfMaxFieldDiff, 0);
    addToRunTimeSelectionTable(functionObject, stopIfMaxFieldDiff, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::stopIfMaxFieldDiff::stopIfMaxFieldDiff
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    field1Ptr_(nullptr),
    field2Ptr_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::stopIfMaxFieldDiff::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {

        field1Name_ = dict.get<word>("field1");
        field2Name_ = dict.get<word>("field2");

        return true;
    }

    return false;
}


bool Foam::functionObjects::stopIfMaxFieldDiff::execute()
{
    return true;
}


bool Foam::functionObjects::stopIfMaxFieldDiff::write()
{   
    //- Set pointers
    if (field1Ptr_ == nullptr)
    {
        field1Ptr_ = &mesh_.lookupObject<volScalarField>(field1Name_);
    }
    if (field2Ptr_ == nullptr)
    {
        field2Ptr_ = &mesh_.lookupObject<volScalarField>(field2Name_);
    }
    const volScalarField& field1(*field1Ptr_);
    const volScalarField& field2(*field2Ptr_);

    //- Check
    if (max(field1-field2).value() > 0.0)
    {
        Log << "    Terminated by " << type() << " " << name() << " at time = "
            << mesh_.time().timeName() << " s" << endl;

        Time& time(const_cast<Time&>(mesh_.time()));
        time.writeAndEnd();
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
