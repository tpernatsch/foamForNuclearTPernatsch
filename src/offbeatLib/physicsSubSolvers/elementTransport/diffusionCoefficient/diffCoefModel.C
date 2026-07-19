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

#include "diffCoefModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diffCoefModel, 0);
    defineRunTimeSelectionTable(diffCoefModel, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffCoefModel::diffCoefModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word defaultModel
)
:
    mesh_(mesh),
    crackFp_(dict.lookupOrDefault<bool>("crackFpRelease", false)),
    failureTime_(dict.lookupOrDefault<scalar>("failureTime", GREAT)),
    FpNamesMat_()
{
  if (dict.found("fissionProductsTransport"))
  {
    const dictionary& FpDict = dict.subDict("fissionProductsTransport");
    FpNamesMat_ = wordList(FpDict.lookup("fissionProductsName"));
    d1_.resize(FpNamesMat_.size());
    q1_.resize(FpNamesMat_.size());
    d2_.resize(FpNamesMat_.size());
    q2_.resize(FpNamesMat_.size());

    forAll(FpNamesMat_, i)
    {
      const dictionary& coefDict = FpDict.subDict(FpNamesMat_[i]);
      d1_[i] = coefDict.lookupOrDefault<scalar>("d1", 0);
      q1_[i] = coefDict.lookupOrDefault<scalar>("q1", 0);
      d2_[i] = coefDict.lookupOrDefault<scalar>("d2", 0);
      q2_[i] = coefDict.lookupOrDefault<scalar>("q2", 0);
    }
  }

}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::diffCoefModel>
Foam::diffCoefModel::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word defaultModel
)
{
    word diffCoefModelName;

    diffCoefModelName =
    dict.lookupOrDefault<word>("diffCoefModel", defaultModel);

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(diffCoefModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("diffCoefModel::New(const fvMesh&, const dictionary&, const word)")
            << "Unknown diffCoefModel type "
            << diffCoefModelName << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting diffCoefModel type "
            << diffCoefModelName << endl;
    }

    return autoPtr<diffCoefModel>(cstrIter()(mesh, dict, defaultModel));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffCoefModel::~diffCoefModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
