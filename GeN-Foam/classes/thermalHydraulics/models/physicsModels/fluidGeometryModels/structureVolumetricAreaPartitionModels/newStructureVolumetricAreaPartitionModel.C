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

#include "structureVolumetricAreaPartitionModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::structureVolumetricAreaPartitionModel> 
Foam::structureVolumetricAreaPartitionModel::New
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const fluid& fluid1,
    const fluid& fluid2,
    const structureModel& structure
)
{
    word type(dict.lookup("type"));

    Info<< "Selecting structureVolumetricAreaPartitionModel : " 
        << type << endl;

    structureVolumetricAreaPartitionModelsConstructorTable::iterator 
        cstrIter = 
        structureVolumetricAreaPartitionModelsConstructorTablePtr_->find
        (
            type
        );

    if 
    (
        cstrIter 
        == 
        structureVolumetricAreaPartitionModelsConstructorTablePtr_->end()
    )
    {
        FatalErrorInFunction
            << "Unknown structureVolumetricAreaPartitionModel type "
            << type << endl << endl
            << "Valid structureVolumetricAreaPartitionModel types are : " 
            << endl
            << structureVolumetricAreaPartitionModelsConstructorTablePtr_
            ->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(objReg, dict, fluid1, fluid2, structure);
}


// ************************************************************************* //
