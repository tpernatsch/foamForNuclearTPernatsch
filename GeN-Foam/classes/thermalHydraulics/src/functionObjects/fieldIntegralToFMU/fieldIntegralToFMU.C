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
#if __has_include("commDataLayer.H") 
#include "commDataLayer.H"
#define isCommDataLayerIncluded
#endif

#ifdef isCommDataLayerIncluded


#include "fieldIntegralToFMU.H"
#include "addToRunTimeSelectionTable.H"
#include "commDataLayer.H"
#include "externalIOObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
/*
// public fvMeshFunctionObject  //This would also work
// but one would have to put the functionObject in
// the controlDict instead of the externalCouplingDict
// which goes agaist the logif of FMU4FOAM
namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldIntegralToFMU, 0);
    addToRunTimeSelectionTable(functionObject, fieldIntegralToFMU, dictionary);
}
}
*/

namespace Foam
{
namespace externalIOObject
{
    defineTypeNameAndDebug(fieldIntegralToFMU, 0);
    addToRunTimeSelectionTable(externalIOObject, fieldIntegralToFMU, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalIOObject::fieldIntegralToFMU::fieldIntegralToFMU
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    externalIOObject(name, runTime, dict),
    fieldName_(dict.get<word>("fieldName")),
    mesh_
    (
        refCast<const fvMesh>
        (
            time_.lookupObject<objectRegistry>
            (
                dict.getOrDefault("region", polyMesh::defaultRegion)
            )
        )
    ),
    cellZone_(dict.get<word>("cellZone")),
    nameFMU_(dict.get<word>("nameFMU")),    
    fieldPtr_(nullptr)
{
    read(dict);
    execute();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::externalIOObject::fieldIntegralToFMU::read(const dictionary& dict)
{

    commDataLayer& data = commDataLayer::New(time_);

    data.storeObj(0.0,nameFMU_,commDataLayer::causality::out);
    
    return false;
}

bool Foam::externalIOObject::fieldIntegralToFMU::execute()
{
    commDataLayer& data = commDataLayer::New(time_);
    
    scalar& result = data.getObj<scalar>(nameFMU_,commDataLayer::causality::out);

    const volScalarField& field = mesh_.lookupObject<volScalarField>(fieldName_);

    label cellZoneID = mesh_.cellZones().findZoneID(cellZone_);
    const cellZone& tgtCellZone = mesh_.cellZones()[cellZoneID];

    scalarField fieldZone(field,tgtCellZone);
    scalarField  volZone(mesh_.V(),tgtCellZone);

    result = gSum(fieldZone * volZone);

    return false;
}

bool Foam::externalIOObject::fieldIntegralToFMU::write()
{ 
    return false;
}

#endif
// ************************************************************************* //
