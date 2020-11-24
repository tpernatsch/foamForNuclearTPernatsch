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

#include "regimeMapModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regimeMapModel, 0);
    defineRunTimeSelectionTable
    (
        regimeMapModel, 
        regimeMapModels
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regimeMapModel::regimeMapModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& physicsModelsDict
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName("regimeMapModel", typeName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),    
    mesh_(mesh),
    physicsModelsDict_(physicsModelsDict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regimeMapModel::setRequiresModelCorrection()
{
    //- I cannot do both these actions in a single loop. Imagine having 3 
    //  regimes, r0,r1,r2, of which r1 is an interpolated regime between
    //  r0 and r2 (its name will be r0.r2 but for the sake of discussiong
    //  let us call it simply r1). Imagine r0 and r1 are the only regimes
    //  present. It means that r0 requires model correction as well as r1.
    //  Since r1 is interpolated from r0 and r2, correcting the models of r1
    //  consists of correcting the models of both r0 and r2. The model
    //  correction is done elsewhere in the code and is done on the basis
    //  of the requiresModelCorrection flag, set here. If both the ifs are
    //  handled in only one loop, if the order in which the regimes are
    //  iterated is r0,r1,r2, it is easy to see that, at the end of the loop,
    //  the requiresModelCorrection of r2 will be set to false, even though
    //  it should be true, as r2 models are used by the interpolated regime r1.
    //  Thus, the loop is split in two, first to set the flag based on
    //  the physical existence of the regime, secondly to account for the
    //  fact that interpolated regimes rely on other models in regimes that
    //  might not be present (yet that need to be updated for interpolation
    //  purposes)
    forAllIter
    (
        regimeTable,
        regimes_,
        iter
    )
    {
        regime& regime(iter()());
        regime.requiresModelCorrection() = regime.isCurrentlyPresent();
    }

    forAllIter
    (
        regimeTable,
        regimes_,
        iter
    )
    {
        regime& regime(iter()());
        if (regime.isInterpolated() and regime.isCurrentlyPresent())
        {
            regime.regime1().requiresModelCorrection() = true;
            regime.regime2().requiresModelCorrection() = true;
        }
    }
}

// ************************************************************************* //
