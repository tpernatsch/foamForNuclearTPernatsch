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

#include "fluidInterfacialAreaModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::fluidInterfacialAreaModel> 
Foam::fluidInterfacialAreaModel::New
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const fluid& fluid1,
    const fluid& fluid2
)
{
    word type(dict.lookup("type"));

    Info<< "Selecting fluidInterfacialAreaModel: " 
        << type << endl;

    fluidInterfacialAreaModelsConstructorTable::iterator cstrIter =
        fluidInterfacialAreaModelsConstructorTablePtr_->find
        (
            type
        );

    if (cstrIter == fluidInterfacialAreaModelsConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fluidInterfacialAreaModelType type "
            << type << endl << endl
            << "Valid fluidInterfacialAreaModel types are : " << endl
            << fluidInterfacialAreaModelsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(objReg, dict, fluid1, fluid2);
}


// ************************************************************************* //
