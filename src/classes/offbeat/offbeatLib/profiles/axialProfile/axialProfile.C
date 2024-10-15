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

#include "axialProfile.H"
#include "zeroGradientFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(axialProfile, 0);
    defineRunTimeSelectionTable(axialProfile, dictionary);

    addToRunTimeSelectionTable
    (
        axialProfile, 
        axialProfile, 
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::axialProfile::axialProfile
(
    const fvMesh& mesh,    
    const dictionary& dict,
    const labelList& addr,
    const vector& axialDirection,
    scalar zMin,
    scalar zMax 
)
:
    mesh_(mesh),
    addr_(addr),
    profile_(addr_.size(), 1.0),
    axialDirection_(axialDirection),
    zMin_(zMin),
    zMax_(zMax)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::axialProfile>
Foam::axialProfile::New
(
    const fvMesh& mesh,    
    const dictionary& dict,
    const labelList& addr,
    const vector& axialDirection,
    scalar zMin,
    scalar zMax 
)
{
    word type = dict.lookupOrDefault<word>("type", "flat");

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {                
        FatalErrorIn("axialProfile::New(const fvMesh&, const dictionary&, const labelList&)")
            << "Unknown axial profile type "
            << type << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }
    
    Info<< "Selecting axial profile " << type << endl;
    
    return autoPtr<axialProfile>(cstrIter()
        (mesh, dict, addr, axialDirection, zMin, zMax));
     
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::axialProfile::~axialProfile()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::axialProfile::axialLocations() const
{    
    const fvMesh& referenceMesh = 
    ( 
        mesh_.foundObject<fvMesh>("referenceMesh")
    ) ?
    mesh_.lookupObject<fvMesh>("referenceMesh")
    : 
    mesh_;
    
    tmp<scalarField> tz(new scalarField(addr_.size()));
    scalarField& z = tz.ref();
    
    const vectorField& C = referenceMesh.C();
    
    forAll(addr_, addrI)
    {
        const label cellI = addr_[addrI];
        z[addrI] = ((C[cellI] & axialDirection_)-zMin_)/(zMax_-zMin_);
    }
    
    return tz;
}


void Foam::axialProfile::checkNormalization() const
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
        WarningIn("Foam::axialProfile::checkNormalization()")
            << "Axial profile is not normalized, with an error of " 
            << error << "%" << nl << endl;
    }
}

// ************************************************************************* //
