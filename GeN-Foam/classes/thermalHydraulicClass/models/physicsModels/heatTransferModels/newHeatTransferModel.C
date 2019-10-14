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

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::heatTransferModel> Foam::heatTransferModel::New
(
    const objectRegistry& objReg,
    const dictionary& heatTransferModelsDict,
    const FFPair& FFPair,
    const word nameBulk,
    const word nameInterface
)
{
    word pair(nameBulk+"."+nameInterface);

    const dictionary& dict(heatTransferModelsDict.subDict(pair));

    word type(dict.lookup("type"));

    Info<< "Selecting heatTransferModel btw. bulk of "
        << nameBulk << " interfacing with " << nameInterface
        <<  ": " << type << endl;

    FFHeatTransferModelsConstructorTable::iterator cstrIter =
        FFHeatTransferModelsConstructorTablePtr_->find(type);

    if (cstrIter == FFHeatTransferModelsConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fluid-fluid heatTransferModel type "
            << type << endl << endl
            << "Valid heatTransferModel types are: " << endl
            << FFHeatTransferModelsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(objReg, dict, FFPair, nameBulk, nameInterface);
}

Foam::autoPtr<Foam::heatTransferModel> Foam::heatTransferModel::New
(
    const objectRegistry& objReg,
    const dictionary& heatTransferModelsDict,
    const FSPair& FSPair,
    const word region
)
{
    word pair;
    word p1p2(FSPair.nameFluid()+"."+region);
    word p2p1(region+"."+FSPair.nameFluid());
    
    pair =  (heatTransferModelsDict.found(p1p2)) 
            ? 
            p1p2 
            : 
            ( 
                (heatTransferModelsDict.found(p2p1)) 
                ? 
                p2p1 
                : 
                word("") 
            );

    if (pair == "")
    {
        FatalErrorInFunction
            << "Fluid-structure heatTransferModel for "
            << FSPair.nameFluid() << " and " << region
            << " not found!" << endl
            << exit(FatalError);
    }

    const dictionary& dict(heatTransferModelsDict.subDict(pair));

    pair = p1p2;    //- By my standardized scheme, name1 should always be the
                    //  fluid, while name2 the name of the structure region

    word type(dict.lookup("type"));

    Info<< "Selecting heatTransferModel btw. bulk of "
        << FSPair.nameFluid() << " interfacing with " << region
        <<  ": " << type << endl;

    FSHeatTransferModelsConstructorTable::iterator cstrIter =
        FSHeatTransferModelsConstructorTablePtr_->find(type);

    if (cstrIter == FSHeatTransferModelsConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown type type "
            << type << endl << endl
            << "Valid heatTransferModel types are : " << endl
            << FSHeatTransferModelsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(objReg, dict, FSPair, region);
}

// ************************************************************************* //
