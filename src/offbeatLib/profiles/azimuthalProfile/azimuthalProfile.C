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

#include "azimuthalProfile.H"
#include "zeroGradientFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(azimuthalProfile, 0);
    defineRunTimeSelectionTable(azimuthalProfile, dictionary);

    addToRunTimeSelectionTable
    (
        azimuthalProfile, 
        azimuthalProfile, 
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::azimuthalProfile::azimuthalProfile
(
    const fvMesh& mesh,    
    const dictionary& dict,
    const labelList& addr,
    const vector& axialDirection,
    const vector& radialDirection
)
:
    mesh_(mesh),
    addr_(addr),
    profile_(addr_.size(), 1.0),
    axialDirection_(axialDirection),
    radialDirection_(radialDirection)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::azimuthalProfile>
Foam::azimuthalProfile::New
(
    const fvMesh& mesh,    
    const dictionary& dict,
    const labelList& addr,
    const vector& axialDirection,
    const vector& radialDirection
)
{
    word type = dict.lookupOrDefault<word>("type", "flat");

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("azimuthalProfile::New(const fvMesh&, const dictionary&, const labelList&)")
            << "Unknown axial profile type "
            << type << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }
    
    Info<< "Selecting azimuthal profile --> " << type << endl;
    
    return autoPtr<azimuthalProfile>(cstrIter()
        (mesh, dict, addr, axialDirection, radialDirection));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::azimuthalProfile::~azimuthalProfile()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::azimuthalProfile::azimuthalLocations() const
{
    const fvMesh& referenceMesh = 
    ( 
        mesh_.foundObject<fvMesh>("referenceMesh")
    ) ?
    mesh_.lookupObject<fvMesh>("referenceMesh")
    : 
    mesh_;
    
    tmp<scalarField> ttheta(new scalarField(addr_.size()));
    scalarField& theta = ttheta.ref();
    
    const vectorField& C = referenceMesh.C();
    
    forAll(addr_, addrI)
    {
        const label cellI = addr_[addrI];

        // Calculate position vector and normalize
        vector r = C[cellI] - (C[cellI] & axialDirection_)*axialDirection_;
        r /= mag(r);

        // Ensure that radialdirection is perpendicular to axialDirection
        // and normalize
        vector radialProjectedDirection = radialDirection_ - (
            radialDirection_ & axialDirection_)*axialDirection_;
        radialProjectedDirection /= mag(radialProjectedDirection);

        // Calculate angle and normalize
        theta[addrI] = acos(
            min(1.0, max(-1.0, (r & radialProjectedDirection))));
        theta[addrI] = theta[addrI] / (2 * Foam::constant::mathematical::pi);
    }
    
    return ttheta;
}


void Foam::azimuthalProfile::checkNormalization() const
{

    const fvMesh& referenceMesh = 
    ( 
        mesh_.foundObject<fvMesh>("referenceMesh")
    ) ?
    mesh_.lookupObject<fvMesh>("referenceMesh")
    : 
    mesh_;
 
    scalarField volumes(referenceMesh.V(), addr_);
    scalar error =
    100*
    (
        mag(gSum(volumes) - gSum(volumes*profile_))/gSum(volumes)
    );

    if(error > 0.01)
    {
        WarningIn("Foam::azimuthalProfile::checkNormalization()")
            << "Azimuthal profile is not normalized, with an error of " 
            << error << "%" << nl << endl;
    }
}

// ************************************************************************* //
