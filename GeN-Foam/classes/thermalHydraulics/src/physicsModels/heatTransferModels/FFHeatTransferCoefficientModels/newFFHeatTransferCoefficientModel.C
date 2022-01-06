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

#include "FFPair.H"
#include "FFHeatTransferCoefficientModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::FFHeatTransferCoefficientModel> 
Foam::FFHeatTransferCoefficientModel::New
(
    const FFPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
{
    word type(dict.get<word>("type"));

    Info<< "Selecting FFHeatTransferCoefficientModel for pair " << pair.name() 
        << " on side " << dict.dictName() << ": " << type << endl;

    auto* ctorPtr = FFHeatTransferCoefficientModelsConstructorTable(type);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown fluid-structure FFHeatTransferCoefficientModel type "
            << type << endl << endl
            << "Valid fluid-structure FFHeatTransferCoefficientModel types "
            << "are: " << endl
            << FFHeatTransferCoefficientModelsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return
        autoPtr<FFHeatTransferCoefficientModel>
        (
            ctorPtr(pair, dict, objReg)
        );

}

// ************************************************************************* //
