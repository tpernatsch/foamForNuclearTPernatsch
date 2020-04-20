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

#include "noRegimeMap.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regimeMapModels
{
    defineTypeNameAndDebug(noRegimeMap, 0);
    addToRunTimeSelectionTable
    (
        regimeMapModel, 
        noRegimeMap, 
        regimeMapModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regimeMapModels::noRegimeMap::noRegimeMap
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& physicsModelsDict
)
:
    regimeMapModel
    (
        mesh,
        dict,
        physicsModelsDict
    )

{
    //- Create a default regime. It requires not special subDict in
    //  physicsModelsDict_
    word regimeName("defaultRegime");
    regimes_.set
    (
        "defaultRegime",
        autoPtr<regime>
        (
            new regime
            (
                mesh_,
                regimeName,
                this->physicsModelsDict_
            )
        )
    );

    regime& defaultRegime(regimes_["defaultRegime"]());
    
    //- Set cellField to all 1s (regime exists everywhere)
    defaultRegime.cellField() = scalarField(mesh_.cells().size(), 1.0);
    
    //- Set cellList to list of all mesh cell indices (regime exists 
    //  everywhere)
    labelList& cellList(defaultRegime.cellList());
    cellList = labelList(0); //- Reset just to be sure
    forAll(mesh_.cells(), i)
    {
        cellList.append(i);
    }

    //- Set this flag to true so that the physicsModels within the regime will
    //  get always corrected
    defaultRegime.requiresModelCorrection() = true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regimeMapModels::noRegimeMap::correct()
{}


// ************************************************************************* //
