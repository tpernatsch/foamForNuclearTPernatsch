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

#include "heatTransferModel.H"
#include "myStringOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatTransferModel, 0);
    defineRunTimeSelectionTable(heatTransferModel, FFHeatTransferModels);
    defineRunTimeSelectionTable(heatTransferModel, FSHeatTransferModels);
}

const Foam::dimensionSet Foam::heatTransferModel::dimHtc(1, 0, -3, 1, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModel::heatTransferModel
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FFPair& FFPair,
    const word nameBulk,
    const word nameInterface
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
    pairName_(nameBulk+"."+nameInterface),
    nameBulk_(nameBulk),
    nameInterface_(nameInterface),
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

Foam::heatTransferModel::heatTransferModel
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FSPair& FSPair,
    const wordList& regions
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName(typeName, FSPair.name()+"."+dict.name()),
            FSPair.mesh().time().timeName(),
            objReg,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(FSPair.mesh()),
    pairName_
    (
        FSPair.name()
    ),
    nameBulk_(FSPair.fluidName()),
    nameInterface_(dict.dictName()),
    FFPair_(nullptr),
    FSPair_(&FSPair),
    cellList_(0),
    cellField_(FSPair.mesh().cells().size(), 0.0)
{    
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
Foam::heatTransferModel::makeInTable
(
    const objectRegistry& srcObjReg,
    const objectRegistry& dstObjReg,
    const word key,
    const dictionary& dict,
    HashTable<autoPtr<heatTransferModel>, word, word::hash>& table
)
{
    wordList entries = myStringOps::split<word>(key, '.');
    word entry1(entries[0]);
    word entry2((entries.size() == 1) ? entries[0] : entries[1]);

    int nFluids(srcObjReg.lookupClass<FSPair>().size());
    
    if (nFluids != 1)
    {
        //- If the key pertains to a fluid-structure model
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
                    const dictionary& subDictFS(dictFS.subDict(key));

                    //- Determine possible multiple regions from dictName. If 
                    //  the key looks something like "region0:region1:huhu"  
                    //  then regions will be = ["region0", "region1", "huhu"].
                    //  Also works for keys that consist of one region only
                    const wordList& regions
                    (
                        myStringOps::split<word>(subDictFS.dictName(), ':')
                    );
                    table.insert
                    (
                        IOobject::groupName(pair.name(), key),
                        heatTransferModel::New
                        (
                            dstObjReg,
                            subDictFS,
                            pair,
                            regions
                        )
                    );
                }
            }
            else
            {
                FatalErrorInFunction
                    << "Key " << key << " contains unkown phases/regions!" << endl
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

                //- Fluid-fluid heatTransferModels are two-resistance models, i.e.
                //  you need two heat transfer coeffs, one from the bulk of each
                //  phase to the phase interface. These are specified by the key
                //  ordering, i.e. the model under the key "entry1.entry2" will
                //  be applied between the bulk of fluid 1 and the interface,
                //  while "entry2.entry1" will be applied from the bulk of fluid 2
                //  to the interface. Such ordering does not matter for fluid-
                //  structure heat transfer models for obvious reasons
                table.insert
                (
                    key12,
                    heatTransferModel::New
                    (
                        dstObjReg,
                        dict.subDict(key12),
                        pair,
                        entry1,
                        entry2 
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

            const dictionary& subDict(dict.subDict(key));

            //- Determine possible multiple regions from dictName. If 
            //  the key looks something like "region0:region1:huhu"  
            //  then regions will be = ["region0", "region1", "huhu"].
            //  Also works for keys that consist of one region only
            const wordList& regions
            (
                myStringOps::split<word>(subDict.dictName(), ':')
            );

            table.insert
            (
                IOobject::groupName(pair.name(), key),
                heatTransferModel::New
                (
                    dstObjReg,
                    subDict,
                    pair,
                    regions
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// ************************************************************************* //
