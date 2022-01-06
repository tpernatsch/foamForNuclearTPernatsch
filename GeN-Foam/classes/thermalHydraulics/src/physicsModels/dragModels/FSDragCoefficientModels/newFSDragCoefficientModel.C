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

#include "FSPair.H"
#include "FSDragCoefficientModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::FSDragCoefficientModel> Foam::FSDragCoefficientModel::New
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
{
    word type(dict.get<word>("type"));

    //- Set fluidName to "fluid" if the fluid has no name (just for Info
    //  cleanliness when using the monoPhase solver)
    word fluidName(pair.fluidName() == "" ? "fluid" : pair.fluidName());

    Info<< "Selecting FSDragCoefficientModel for pair " << pair.name()
        << ": " << type << endl;

    auto* ctorPtr = FSDragCoefficientModelsConstructorTable(type);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown fluid-structure FSDragCoefficientModel type "
            << type << endl << endl
            << "Valid fluid-structure FSDragCoefficientModel types are: " << endl
            << FSDragCoefficientModelsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return
        autoPtr<FSDragCoefficientModel>
        (
            ctorPtr(pair, dict, objReg)
        );
}

// ************************************************************************* //
