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
#include "dispersionModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dispersionModel> Foam::dispersionModel::New
(
    const FFPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
{
    word type(dict.lookup("type"));

    Info<< "Selecting dispersionModel for pair " << pair.name() << ": " 
        << type << endl;

    dispersionModelsConstructorTable::iterator 
        cstrIter = dispersionModelsConstructorTablePtr_->find(type);

    if 
    (
        cstrIter 
        == 
        dispersionModelsConstructorTablePtr_->end()
    )
    {
        FatalErrorInFunction
            << "Unknown dispersionModel type "
            << type << endl << endl
            << "Valid dispersionModel types are : " 
            << endl
            << dispersionModelsConstructorTablePtr_
            ->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(pair, dict, objReg);
}


// ************************************************************************* //
