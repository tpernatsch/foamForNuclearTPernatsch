/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "constantBurnup.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchField.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantBurnup, 0);
    addToRunTimeSelectionTable
    (
        burnup,
        constantBurnup,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantBurnup::constantBurnup
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& burnupDict
)
:
    burnup(mesh, mat, burnupDict),
    Bu_
    (
        IOobject
        (
            "Bu",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Bu", dimless, 0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    maxBuIncrease_(GREAT)
{
    // Prepare old time field
    Bu_.oldTime();

    // Read maximum allowed burnup increment from controlDict
    const dictionary& controlDict(mesh_.time().controlDict());

    bool adjustTime = controlDict.lookupOrDefault<bool>("adjustableTimeStep", false);

    if(adjustTime)
    {
        maxBuIncrease_ = controlDict.lookupOrDefault<scalar>("maxBurnupIncrease", GREAT);

        // Transform from MWd/kg to MWd/t (Bu unit in OFFBEAT)
        maxBuIncrease_ *= 1000;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantBurnup::~constantBurnup()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::constantBurnup::nextDeltaT()
{
    scalar nextDeltaT(0);
    scalar currentDeltaT = mesh_.time().deltaT().value();

    const scalarField& oldBui(Bu_.oldTime().internalField());
    const scalarField& currentBui(Bu_.internalField());

    // Calculate deltaT multiplier based on maximum allowed burnup increase
    scalar burnupIncrement(gMax(currentBui - oldBui));
    scalar deltaTMultiplier = maxBuIncrease_ / (burnupIncrement + SMALL);

    // Set deltaT
    nextDeltaT = deltaTMultiplier*currentDeltaT;

    Info<< "Maximum deltaT calculated by burnup model: "
        << mesh_.time().timeToUserTime(nextDeltaT)
        << " with a maximum burnup increment of " << burnupIncrement
        << endl;

    return nextDeltaT;
}
// ************************************************************************* //
