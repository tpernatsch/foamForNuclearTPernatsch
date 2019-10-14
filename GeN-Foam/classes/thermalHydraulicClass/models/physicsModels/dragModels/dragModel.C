/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "dragModel.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dragModel, 0);
    defineRunTimeSelectionTable(dragModel, FFDragModels);
    defineRunTimeSelectionTable(dragModel, FSDragModels);
}

const Foam::dimensionSet Foam::dragModel::dimK(1, -3, -1, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModel::dragModel
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FFPair& FFPair
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName(typeName, FFPair.name()),
            FFPair.mesh().time().timeName(),
            objReg,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(FFPair.mesh()),
    pairName_(FFPair.name()),
    name1_(FFPair.name1()),
    name2_(FFPair.name2()),
    withStructure_(false)
{}

Foam::dragModel::dragModel
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FSPair& FSPair,
    const word region
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName(typeName, FSPair.name1()+"."+region),
            FSPair.mesh().time().timeName(),
            objReg,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(FSPair.mesh()),
    pairName_(FSPair.name1()+"."+region),
    name1_(FSPair.name1()),
    name2_(region),
    withStructure_(true)
{}

// * * * * * * * * * * * * * * * Static Member Functions   * * * * * * * * * //

void
Foam::dragModel::makeInTable
(
    const objectRegistry& srcObjReg,
    const objectRegistry& dstObjReg,
    const word key,
    const dictionary& dict,
    HashTable<autoPtr<dragModel>, word, word::hash>& table
)
{
    word phase1Name(word(stringOps::split(key, ".").first()));
    word phase2Name(word(stringOps::split(key, ".").last()));
    word key12(key);
    word key21(phase2Name+"."+phase1Name);
    word key1s(phase1Name+".structure");
    word key2s(phase2Name+".structure");
    bool fatalError(false);
    
    if 
    (
        srcObjReg.foundObject<FFPair>(key12) 
        or 
        srcObjReg.foundObject<FFPair>(key21)
    )
    {
        const FFPair& pair
        (
            (srcObjReg.names().found(key12)) ?
            srcObjReg.lookupObject<FFPair>(key12) :
            srcObjReg.lookupObject<FFPair>(key21)
        );
        table.insert
        (
            pair.name(),
            dragModel::New
            (
                dstObjReg,
                dict,
                pair 
            )
        );
    }
    else if 
    (
        srcObjReg.foundObject<FSPair>(key1s) 
        or 
        srcObjReg.foundObject<FSPair>(key2s)
    )
    {
        word keyStruct;
        word keyRegion;
        word region;
        if (srcObjReg.names().found(key1s))
        {   
            keyStruct = key1s;
            keyRegion = key12;
            region = phase2Name;
        }
        else
        {
            keyStruct = key2s;
            keyRegion = key21;
            region = phase1Name;
        }
        const FSPair& pair
        (
            srcObjReg.lookupObject<FSPair>(keyStruct)
        );
        if (!pair.structureRegions().found(region))
        {
            fatalError = true;
        }
        table.insert
        (
            keyRegion,
            dragModel::New
            (
                dstObjReg,
                dict,
                pair,
                region
            )
        );
    }
    else
    {
        fatalError = true;
    }
    if (fatalError)
    {
        FatalErrorInFunction
            << "Key " << key << " contains unkown phases/regions!" << endl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// ************************************************************************* //
