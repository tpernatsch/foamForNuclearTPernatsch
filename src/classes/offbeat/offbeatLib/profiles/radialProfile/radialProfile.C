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

#include "radialProfile.H"
#include "zeroGradientFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(radialProfile, 0);
    defineRunTimeSelectionTable(radialProfile, dictionary);

    addToRunTimeSelectionTable
    (
        radialProfile, 
        radialProfile, 
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radialProfile::radialProfile
(
    const fvMesh& mesh,    
    const dictionary& dict,
    const labelList& addr,
    const vector& axialDirection,
    scalar rMax 
)
:
    mesh_(mesh),
    addr_(addr),
    profile_(addr_.size(), 1.0),
    axialDirection_(axialDirection),
    rMax_(rMax)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::radialProfile>
Foam::radialProfile::New
(
    const fvMesh& mesh,    
    const dictionary& dict,
    const labelList& addr,
    const vector& axialDirection,
    scalar rMax
)
{
    word type = dict.lookupOrDefault<word>("type", "flat");

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("radialProfile::New(const fvMesh&, const dictionary&, const labelList&)")
            << "Unknown radial profile type "
            << type << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }
    
    Info<< "Selecting radial profile --> " << type << endl;
    
    return autoPtr<radialProfile>(cstrIter()(mesh, dict, addr, axialDirection, rMax));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radialProfile::~radialProfile()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::radialProfile::radialLocations() const
{    
    const fvMesh& referenceMesh = 
    ( 
        mesh_.foundObject<fvMesh>("referenceMesh")
    ) ?
    mesh_.lookupObject<fvMesh>("referenceMesh")
    : 
    mesh_;

    tmp<scalarField> tr(new scalarField(addr_.size()));
    scalarField& r = tr.ref();
    
    const vectorField& C = referenceMesh.C();
    
    forAll(addr_, addrI)
    {
        const label cellI = addr_[addrI];
        const vector& p = C[cellI];
        vector n = p - (p & axialDirection_)*axialDirection_;
        r[addrI] = mag(n)/rMax_;
    }
    
    return tr;
}


void Foam::radialProfile::checkNormalization() const
{
    scalarField volumes(mesh_.V(), addr_);
    scalar error =
    100*
    (
        mag(gSum(volumes) - gSum(volumes*profile_))/gSum(volumes)
    );

    if(error > 0.01)
    {
        WarningIn("Foam::radialProfile::checkNormalization()")
            << "Radial profile is not normalized, with an error of " 
            << error << "%" << endl;
    }
}

// ************************************************************************* //
