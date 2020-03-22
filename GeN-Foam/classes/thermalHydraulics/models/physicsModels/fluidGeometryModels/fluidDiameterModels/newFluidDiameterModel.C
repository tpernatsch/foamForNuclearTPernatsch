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

#include "fluidDiameterModel.H"
#include "fluid.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidDiameterModel> Foam::fluidDiameterModel::New
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const fluid& fluid
)
{
    word type(dict.lookup("type"));

    Info<< "Selecting fluidDiameterModel for dispersed fluid "
        << fluid.name() << ": " << type << endl;

    fluidDiameterModelsConstructorTable::iterator cstrIter =
        fluidDiameterModelsConstructorTablePtr_->find(type);

    if (cstrIter == fluidDiameterModelsConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown type type "
            << type << endl << endl
            << "Valid fluidDiameterModel types are : " << endl
            << fluidDiameterModelsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(objReg, dict, fluid);
}


// ************************************************************************* //
