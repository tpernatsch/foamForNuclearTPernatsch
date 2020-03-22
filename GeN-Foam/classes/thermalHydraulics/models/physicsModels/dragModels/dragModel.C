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
#include "myStringOps.H"

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
            IOobject::groupName("dragModel."+typeName, FFPair.name()),
            FFPair.mesh().time().timeName(),
            objReg,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(FFPair.mesh()),
    pairName_(FFPair.name()),
    FFPair_(&FFPair),
    FSPair_(nullptr),
    cellList_(0),
    cellField_(FFPair.mesh().cells().size(), 0.0)
{
    forAll(mesh_.cells(), i)
    {
        cellList_.append(i);
        cellField_[i] = 1.0;
    }
}

Foam::dragModel::dragModel
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FSPair& FSPair
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName("dragModel."+typeName, FSPair.name()),
            FSPair.mesh().time().timeName(),
            objReg,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(FSPair.mesh()),
    pairName_(FSPair.name()),
    FFPair_(nullptr),
    FSPair_(&FSPair),
    cellList_(0),
    cellField_(FSPair.mesh().cells().size(), 0.0)
{
    //- Determine possible multiple regions from dictName. If the key looks
    //  something like "region0:region1:huhu" then regions will be 
    //  = ["region0", "region1", "huhu"]. Also works for keys that consist of
    //  one region only
    wordList regions(myStringOps::split<word>(dict.dictName(), ':'));
    
    forAll(regions, i)
    {
        word region(regions[i]);
        const labelList& cellList(FSPair_->structure().cellLists()[region]);
        forAll(cellList, j)
        {
            label cellj(cellList[j]);
            cellList_.append(cellj);
            cellField_[cellj] = 1.0;
        }
    }
}


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
    wordList entries = myStringOps::split<word>(key, '.');
    word entry1(entries[0]);
    word entry2((entries.size() == 1) ? entries[0] : entries[1]);

    int nFluids(srcObjReg.lookupClass<fluid>().size());
    
    if (nFluids != 1)
    {
        if (entry1 == "structure" or entry2 == "structure")
        {
            word keyFS = 
                (entry1 == "structure") ?
                IOobject::groupName(entry2, entry1) :
                IOobject::groupName(entry1, entry2);
            
            if
            (
                srcObjReg.foundObject<FSPair>(keyFS) 
            )
            {
                const FSPair& pair
                (
                    srcObjReg.lookupObject<FSPair>(keyFS)
                );
                const dictionary& dictFS(dict.subDict(key));
                
                forAllConstIter
                (
                    dictionary,
                    dictFS,
                    iter
                )
                {
                    word key((*iter).keyword());
                    table.insert
                    (
                        IOobject::groupName(pair.name(), key),
                        dragModel::New
                        (
                            dstObjReg,
                            dictFS.subDict(key),
                            pair
                        )
                    );
                }
            }
            else
            {
                FatalErrorInFunction
                    << "Key " << key << " contains unkown phases/regions!" 
                    << exit(FatalError);
            }
        }
        else //- key is of a fluid-fluid pair, much simpler to handle
        {
            word key12(IOobject::groupName(entry1, entry2));
            word key21(IOobject::groupName(entry2, entry1));
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
            else
            {
                FatalErrorInFunction
                << "Key " << key << " contains unkown phases/regions!" << endl
                << exit(FatalError);
            } 
        }
    }
    else
    {
        word keyFS = "FSPair";
        if
        (
            srcObjReg.foundObject<FSPair>(keyFS) 
        )
        {
            const FSPair& pair
            (
                srcObjReg.lookupObject<FSPair>(keyFS)
            );
            
            table.insert
            (
                IOobject::groupName(pair.name(), key),
                dragModel::New
                (
                    dstObjReg,
                    dict.subDict(key),
                    pair
                )
            );
        }
        else
        {
            FatalErrorInFunction
                << "Key " << key << " contains unkown phases/regions!" 
                << exit(FatalError);
        }
    }
    
    /*
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
    //- The following only happens for the monoPhase solver
    //  (see constructor of FSPair in FSPair.C)
    else if (srcObjReg.foundObject<FSPair>("FSPair"))
    {
        const FSPair& pair
        (
            srcObjReg.lookupObject<FSPair>("FSPair")
        );
        if (!pair.structureRegions().found(key))
        {
            fatalError = true;
        }
        table.insert
        (
            key,
            dragModel::New
            (
                dstObjReg,
                dict,
                pair,
                key
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
    */
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// ************************************************************************* //
