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

#include "HDiffusionCoefficientModel.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"
#include "materials.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HDiffusionCoefficientModel, 0);
    defineRunTimeSelectionTable(HDiffusionCoefficientModel, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
void Foam::HDiffusionCoefficientModel::correctEffectiveDiffusionCoefficient
(
    const labelList& addr
)
{

    const volScalarField& hydrideVolFrac_(mesh_.lookupObject<volScalarField>("hydridevolumeFraction"));

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        diffH_[cellI] *= 1 - (hydrideVolFrac_[cellI] >= 0.99 ? 1 : hydrideVolFrac_[cellI]);

        const cell& c = mesh_.cells()[cellI];  
        forAll(c, faceI)
        {
            const label patchID = 
            mesh_.boundaryMesh().whichPatch(c[faceI]);

            if (patchID > -1 and diffH_.boundaryField()[patchID].size())
            {
                const label faceID = 
                mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                const scalarField& hydrideVolFracP(hydrideVolFrac_.boundaryField()[patchID]);
                scalarField& diffHP(diffH_.boundaryFieldRef()[patchID]);

                diffHP[faceID] *= 1 - (hydrideVolFracP[faceID] >= 0.99 ? 1 : hydrideVolFracP[faceID]);
            }
        }
    }
    diffH_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HDiffusionCoefficientModel::HDiffusionCoefficientModel
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& lawDict
)
:
    mesh_(mesh),
    mat_(mat),
    lawDict_(lawDict),
    diffH_(createOrLookup<scalar>(mesh, "DH", dimless))
{
   
}

// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::HDiffusionCoefficientModel>
Foam::HDiffusionCoefficientModel::New
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& lawDict
)
{
    word HDiffusionCoefficientModelName;
    HDiffusionCoefficientModelName = lawDict.lookupOrDefault<word>("HDiffusionCoefficientModel", "HDiffusionCoefficientBison");

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(HDiffusionCoefficientModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("HDiffusionCoefficientModel::New(const fvMesh&, const dictionary&, const labelList&)")
            << "Unknown HDiffusionCoefficientModel type "
            << HDiffusionCoefficientModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    Info<< tab << "Selecting HDiffusionCoefficientModel type:  "
        << HDiffusionCoefficientModelName << endl;

    return autoPtr<HDiffusionCoefficientModel>(cstrIter()(mesh, mat, lawDict));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::HDiffusionCoefficientModel::~HDiffusionCoefficientModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //


