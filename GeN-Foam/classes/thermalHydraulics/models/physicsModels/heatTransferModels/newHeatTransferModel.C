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
    const dictionary& dict,
    const FFPair& FFPair,
    const word nameBulk,
    const word nameInterface
)
{
    word pair(nameBulk+"."+nameInterface);

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
    const dictionary& dict,
    const FSPair& FSPair
)
{
    word type(dict.get<word>("type"));

    //- Set fluidName to "fluid" if the fluid has no name (just for Info
    //  cleanliness when using the monoPhase solver)
    word fluidName(FSPair.fluidName() == "" ? "fluid" : FSPair.fluidName());

    Info<< "Selecting heatTransferModel btw. bulk of "
        << fluidName << " interfacing with region(s) " << dict.dictName()
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

    return cstrIter()(objReg, dict, FSPair);
}

// ************************************************************************* //
