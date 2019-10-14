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

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dragModel> Foam::dragModel::New
(
    const objectRegistry& objReg,
    const dictionary& dragModelsDict,
    const FFPair& FFPair
)
{
    word pair;
    word name1(FFPair.name1());
    word name2(FFPair.name2());
    word p1p2(name1+"."+name2);
    word p2p1(name2+"."+name1);
    
    pair =  (dragModelsDict.found(p1p2)) 
            ? 
            p1p2 
            : 
            ( 
                (dragModelsDict.found(p2p1)) 
                ? 
                p2p1 
                : 
                word("") 
            );

    if (pair == "")
    {
        FatalErrorInFunction
            << "Fluid-fluid drag model for "
            << name1 << " and " << name2
            << " not found!" << endl
            << exit(FatalError);
    }

    const dictionary& dict(dragModelsDict.subDict(pair));

    word type(dict.lookup("type"));

    Info<< "Selecting dragModel between " 
        << name1 << " and " << name2 << ": " << type << endl;

    FFDragModelsConstructorTable::iterator cstrIter =
        FFDragModelsConstructorTablePtr_->find(type);

    if (cstrIter == FFDragModelsConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fluid-fluid dragModel type "
            << type << endl << endl
            << "Valid fluid-fluid dragModel types are: " << endl
            << FFDragModelsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(objReg, dict, FFPair);
}

Foam::autoPtr<Foam::dragModel> Foam::dragModel::New
(
    const objectRegistry& objReg,
    const dictionary& dragModelsDict,
    const FSPair& FSPair,
    const word region
)
{
    word pair;
    word name1(FSPair.name1());
    word name2(region);
    word p1p2(name1+"."+name2);
    word p2p1(name2+"."+name1);
    
    pair =  (dragModelsDict.found(p1p2)) 
            ? 
            p1p2 
            : 
            ( 
                (dragModelsDict.found(p2p1)) 
                ? 
                p2p1 
                : 
                word("") 
            );

    if (pair == "")
    {
        FatalErrorInFunction
            << "Fluid-structure drag model for "
            << name1 << " and " << name2
            << " not found!" << endl
            << exit(FatalError);
    }

    const dictionary& dict(dragModelsDict.subDict(pair));

    word type(dict.lookup("type"));

    Info<< "Selecting dragModel between " 
        << name1 << " and " << name2 << ": " << type << endl;

    FSDragModelsConstructorTable::iterator cstrIter =
        FSDragModelsConstructorTablePtr_->find(type);

    if (cstrIter == FSDragModelsConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fluid-structure dragModel type "
            << type << endl << endl
            << "Valid fluid-structure dragModel types are : " << endl
            << FSDragModelsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(objReg, dict, FSPair, region);
}

// ************************************************************************* //
