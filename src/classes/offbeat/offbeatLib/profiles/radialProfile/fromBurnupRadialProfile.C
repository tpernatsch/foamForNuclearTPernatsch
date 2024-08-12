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

#include "fromBurnupRadialProfile.H"
#include "addToRunTimeSelectionTable.H"
#include "InterpolateTables.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fromBurnupRadialProfile, 0);

    addToRunTimeSelectionTable
    (
        radialProfile, 
        fromBurnupRadialProfile, 
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fromBurnupRadialProfile::fromBurnupRadialProfile
(
    const fvMesh& mesh,    
    const dictionary& dict,
    const labelList& addr,
    const vector& axialDirection,
    scalar rMax 
)
:
    radialProfile(mesh, dict, addr, axialDirection, rMax)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fromBurnupRadialProfile::~fromBurnupRadialProfile()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fromBurnupRadialProfile::correct()
{
    if(!mesh_.foundObject<volScalarField>("formFactor"))
    {
        FatalErrorIn("Foam::fromBurnupRadialProfile::correct()")
        << "Power form factor not found. Please enable a neutronics "
        << "subsolver to use this feature. " << exit(FatalError);
    }

	scalarField formFactor
	(
		mesh_.lookupObject<volScalarField>("formFactor").internalField(),
		addr_
	);

	forAll(addr_, i)
    {
        profile_[i] = formFactor[i];
    }
}
// ************************************************************************* //
