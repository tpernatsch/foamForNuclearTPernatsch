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
#include "myOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dragModel, 0);
    defineRunTimeSelectionTable(dragModel, FFDragModels);
    defineRunTimeSelectionTable(dragModel, FSDragModels);
}

const Foam::dimensionSet Foam::dragModel::dimK(1, -3, -1, 0, 0);

const Foam::Enum
<
    Foam::dragModel::transverseDragModel
>
Foam::dragModel::transverseDragModelNames_
(
    {
        { 
            transverseDragModel::isotropic, 
            "isotropic" 
        },
        { 
            transverseDragModel::same, 
            "same" 
        },
        { 
            transverseDragModel::Blasius, 
            "Blasius" 
        },
        { 
            transverseDragModel::GunterShaw, 
            "GunterShaw" 
        }/*,
        { 
            transverseDragModel::anotherTransverseDragModel, 
            "anotherTransverseDragModelName" 
        }*/
    }
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct fluid-fluid drag model
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
    cellField_(FFPair.mesh().cells().size(), 0.0),
    transverseDragModel_    //- All the following initializations are useless 
    (                       //  for fluid-fluid drag, they are here just for
        transverseDragModelNames_.get // safety
        (
            "isotropic"
        )
    ),
    isotropic_(true),
    dta0_(0),
    dta1_(1),
    dpa_(2),
    halfAlphaRhoMagU_(nullptr)
{
    forAll(mesh_.cells(), i)
    {
        cellList_.append(i);
        cellField_[i] = 1.0;
    }
}

//- Construct fluid-structure drag model
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
    cellField_(FSPair.mesh().cells().size(), 0.0),
    transverseDragModel_
    (
        transverseDragModelNames_.get
        (
            dict.lookupOrDefault<word>
            (
                "transverseDragModel", 
                "isotropic"
            )
        )
    ),
    isotropic_(true),
    dta0_(0),
    dta1_(1),
    dpa_(2),
    halfAlphaRhoMagU_(nullptr)
{
    //- Determine possible multiple regions from dictName. If the key looks
    //  something like "region0:region1:huhu" then regions will be 
    //  = ["region0", "region1", "huhu"]. Also works for keys that consist of
    //  one region only
    wordList regions(myOps::split<word>(dict.dictName(), ':'));
    
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

    if (transverseDragModel_ != transverseDragModel::isotropic) 
        isotropic_ = false;

    //- Checks
    wordList transverseDragModelNamesNoIso(0, "");
    forAll(transverseDragModelNames_.names(), i)
    {
        word name(transverseDragModelNames_.names()[i]);
        if (name != "isotropic")
            transverseDragModelNamesNoIso.append(name);
    }

    if (this->found("principalAxis") and isotropic_)
    {
        FatalErrorInFunction
            << "A principalAxis was specified but either no "
            << "transverseDragModel or an isotropic transverseDragModel were " 
            << "specified. Available transverseDragModels are: " 
            << transverseDragModelNamesNoIso
            << exit(FatalError);
    }

    if (!this->found("principalAxis") and !isotropic_)
    {
        FatalErrorInFunction
            << "A transverseDragModel different than isotropic was specified "
            << "but no principalAxis was specified. Possible axes keywords "
            << "are localX, localY or localZ. They refer to the local X, Y or "
            << "Z axes of the structure in the cellZones where this dragModel "
            << "is being specified."
            << exit(FatalError);
    }

    //- Determine principal axis, if any
    if (!isotropic_)
    {
        word principalAxis(this->get<word>("principalAxis"));
        if (principalAxis == "localX")
        {
            dta0_ = 1;
            dta1_ = 2;
            dpa_= 0;  
        }
        else if (principalAxis == "localY")
        {
            dta0_ = 2;
            dta1_ = 0;
            dpa_= 1;  
        }
        else if (principalAxis == "localZ")
        {
            dta0_ = 0;
            dta1_ = 1;
            dpa_= 2; 
        }
    }

    //- Init halfAlphaRhoMagU_
    halfAlphaRhoMagU_ = new scalarField(cellList_.size(), 0);
}


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * //

//- Expand this function as you add more transverse drag models
void Foam::dragModel::correctTransverseKd(volTensorField& Kd) const
{
    switch (transverseDragModel_)
    {
        default :
            break;
        case transverseDragModel::Blasius :
            {
                #include "calcTransverseDragBlasius.H"
            }
            break;
        case transverseDragModel::GunterShaw :
            {
                #include "calcTransverseDragGunterShaw.H"
            }
            break;
    }
}


// * * * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * //

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
    wordList entries = myOps::split<word>(key, '.');
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
}


// ************************************************************************* //
