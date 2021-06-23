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

#include "contactPartitionModel.H"
#include "FSPair.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::contactPartitionModel> 
Foam::contactPartitionModel::New
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
{
    word type(dict.lookup("type"));

    Info<< "Selecting contactPartitionModel for pair " << pair.name() << ": " 
        << type << endl;

    contactPartitionModelsConstructorTable::iterator 
        cstrIter = contactPartitionModelsConstructorTablePtr_->find(type);

    if 
    (
        cstrIter 
        == 
        contactPartitionModelsConstructorTablePtr_->end()
    )
    {
        FatalErrorInFunction
            << "Unknown contactPartitionModel type "
            << type << endl << endl
            << "Valid contactPartitionModel types are : " 
            << endl
            << contactPartitionModelsConstructorTablePtr_
            ->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(pair, dict, objReg);
}


// ************************************************************************* //
