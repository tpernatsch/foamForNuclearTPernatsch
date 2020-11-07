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

#include "poolBoilingModel.H"
#include "fluid.H"
#include "structureModel.H"
#include "FSPair.H"
#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr
<
    Foam::heatTransferModels::nucleateBoilingSubModels::poolBoilingModel
> 
Foam::heatTransferModels::nucleateBoilingSubModels::poolBoilingModel::New
(
    const heatTransferModel& htm,
    const objectRegistry& objReg,
    const FSPair& FSPair
)
{
    word type(htm.subDict("poolBoilingModel").get<word>("type"));
    Info<< "Constructing poolBoilingModel of type " << type 
        << " in region(s): " << FSPair.structure().regions() << endl;

    poolBoilingModelsConstructorTable::iterator cstrIter =
        poolBoilingModelsConstructorTablePtr_->find
        (
            type
        );

    if (cstrIter == poolBoilingModelsConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown type type "
            << type << endl << endl
            << "Valid poolBoilingModel types are : " << endl
            << poolBoilingModelsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(htm, objReg, FSPair);
}


// ************************************************************************* //
