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

#include "constitutiveLaw.H"
#include "addToRunTimeSelectionTable.H"

#include "fvm.H"
#include "fvc.H"
#include "calculatedFvPatchFields.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constitutiveLaw, 0);
    defineRunTimeSelectionTable(constitutiveLaw, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveLaw::constitutiveLaw
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    mesh_(mesh),
    lawDict_(lawDict),    
    damageModel_(damageModel::New(mesh, lawDict))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constitutiveLaw::~constitutiveLaw()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::constitutiveLaw::updateStress
(
    volSymmTensorField& sigma, 
    volScalarField& sigmaHyd, 
    volSymmTensorField& sigmaDev,
    const labelList& addr
)
{
    damageModel_->correct(addr);
    const volScalarField& d(damageModel_->d());
    bool applyDamage(!isNull(d));

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];

        // Apply damage only when applicable
        if(applyDamage)
        {
            sigmaHyd[cellI] *= (1 - d[cellI]);
            sigmaDev[cellI] *= (1 - d[cellI]);
        }

        sigma[cellI] = sigmaHyd[cellI]*I + sigmaDev[cellI];

        const cell& c = mesh_.cells()[cellI];         

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);   

            if (patchID > -1 and sigma.boundaryField()[patchID].size())
            {
                const label faceID = 
                mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                scalarField& sigmaHydF = 
                sigmaHyd.boundaryFieldRef()[patchID];

                symmTensorField& sigmaDevF = 
                sigmaDev.boundaryFieldRef()[patchID];

                // Apply damage only when applicable
                if(applyDamage)
                {
                    const scalarField& dF = 
                    d.boundaryField()[patchID];

                    sigmaHydF[faceID] *= (1 - dF[faceID]);
                    sigmaDevF[faceID] *= (1 - dF[faceID]);
                }

                symmTensorField& sigmaF = sigma.boundaryFieldRef()[patchID];

                sigmaF[faceID] = sigmaHydF[faceID]*I + sigmaDevF[faceID];  
            }
        }
    }
}


Foam::autoPtr<Foam::constitutiveLaw>
Foam::constitutiveLaw::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    // Initialize lawName
    word lawName;

    dictionary lawDict
    (
        dict.subOrEmptyDict("rheologyModelOptions")
    );

    dict.lookup("rheologyModel") >> lawName;

    Info << "     Selecting rheology model " << lawName << endl;

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(lawName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("constitutiveLaw::New(const fvMesh&, const dictionary&, const labelList&)")
            << "Unknown constitutiveLaw type "
            << lawName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    Info<< "     Selecting constitutiveLaw type --> "
        << lawName << endl;
        
    if (debug)
    {
        Info<< "Selecting constitutiveLaw type "
            << lawName << endl;
    }

    return autoPtr<constitutiveLaw>(cstrIter()(mesh, lawDict));
}


// ************************************************************************* //


