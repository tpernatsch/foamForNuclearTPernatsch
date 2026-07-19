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

#include "diffCoefArrhenius.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diffCoefArrhenius, 0);
    addToRunTimeSelectionTable
    (
        diffCoefModel,
        diffCoefArrhenius,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffCoefArrhenius::diffCoefArrhenius
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word defaultModel
)
:
    diffCoefModel(mesh, dict, defaultModel)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffCoefArrhenius::~diffCoefArrhenius()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::diffCoefArrhenius::updateCoef
(
    volScalarField& sf,
    const volScalarField& T,
    const labelList& addr,
    const word FpName
)
{
  bool findFp(false);
  // label FpI(-1);
  //If the material is considered broken, then the diffCoefs for all the species
  //are set to 1e-6. //For Iaea benchmark Cases;
  if (isCracked())
  {
    forAll(addr, i)
    {
      const label cellI = addr[i];
      scalar nominalValue;
      nominalValue = 1e-6;
      sf[cellI] = nominalValue;

      const cell& c = mesh_.cells()[cellI];
      forAll(c, faceI)
      {
        const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

        if (patchID > -1 and sf.boundaryField()[patchID].size())
        {
          // Take references to patch fields
          const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
          scalarField& diffP(sf.boundaryFieldRef()[patchID]);
          scalar nominalValueP;
          nominalValueP = 1e-6;
          diffP[faceID] = nominalValueP;
          }
       }
    }
    return;
  }

  forAll(FpNamesMat_, j)
  {
    if (FpNamesMat_[j] == FpName)
    {
      findFp = true;
      forAll(addr, i)
      {
        const label cellI = addr[i];
        scalar nominalValue;
        nominalValue = d1_[j] * exp(-q1_[j]/8.31446/T[cellI]) + d2_[j] * exp(-q2_[j]/8.31446/T[cellI]);
        sf[cellI] = nominalValue;

        const cell& c = mesh_.cells()[cellI];
        forAll(c, faceI)
        {
          const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

          if (patchID > -1 and sf.boundaryField()[patchID].size())
          {
            // Take references to patch fields
            const scalarField& Tp(T.boundaryField()[patchID]);
            const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
            scalarField& diffP(sf.boundaryFieldRef()[patchID]);
            scalar nominalValueP;
            nominalValueP = d1_[j] * exp(-q1_[j]/8.31446/Tp[faceID])
                          + d2_[j] * exp(-q2_[j]/8.31446/Tp[faceID]);
            diffP[faceID] = nominalValueP;
            }
         }
      }
      return;
    }
  }

  if (! findFp)
  {
    FatalErrorIn("Foam::diffCoefArrhenius::updateCoef")
    << "The diffusion coefficient of " << FpName << " is not provided"
    << exit(FatalError);
  }
}





// ************************************************************************* //
