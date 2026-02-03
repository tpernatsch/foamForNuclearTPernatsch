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

#include "diffCoefSiC.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diffCoefSiC, 0);
    addToRunTimeSelectionTable
    (
        diffCoefModel,
        diffCoefSiC,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffCoefSiC::diffCoefSiC
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word defaultModel
)
:
    diffCoefModel(mesh, dict, defaultModel),
    phi_(mesh_.lookupObject<volScalarField>("fastFluence"))
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffCoefSiC::~diffCoefSiC()
{}

void Foam::diffCoefSiC::updateCsCoef
(
  volScalarField& sf,
  const volScalarField& T,
  const labelList& addr
)
{
  const scalarField& phii(phi_.internalField()); // E > 0.18 MeV
  forAll(addr, i)
  {
      const label cellI = addr[i];
      scalar nominalValue;
      nominalValue = 5.5e-14*exp(phii[cellI]*1e4/1e25*1.1/5.0) * exp(-125e3/8.31446/T[cellI])
                   + 1.6e-2 * exp(-514e3/8.31446/T[cellI]);
      sf[cellI] = nominalValue ;

      const cell& c = mesh_.cells()[cellI];
      forAll(c, faceI)
      {
        const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

        if (patchID > -1 and sf.boundaryField()[patchID].size())
        {
          const scalarField& Tp(T.boundaryField()[patchID]);
          const scalarField& phip(phi_.boundaryField()[patchID]);
          const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
          scalarField& diffP(sf.boundaryFieldRef()[patchID]);
          scalar nominalValueP;
          nominalValueP = 5.5e-14*exp(phip[faceID]*1e4/1e25*1.1/5.0) * exp(-125e3/8.31446/Tp[faceID])
                        + 1.6e-2 * exp(-514e3/8.31446/Tp[faceID]);
          diffP[faceID] = nominalValueP;
        }
      }
    }
}

void Foam::diffCoefSiC::updateKrCoef
(
  volScalarField& sf,
  const volScalarField& T,
  const labelList& addr
)
{
  forAll(addr, i)
  {
      const label cellI = addr[i];
      scalar nominalValue;
      if (T[cellI] < 1673.15)//Parfume:1626K, INL report: 1673.15K
      {
        nominalValue = 8.6e-10 * exp(-326e3/8.31446/T[cellI]);
      }
      else
      {
        nominalValue = 3.7e1 * exp(-657e3/8.31446/T[cellI]);
      }
      sf[cellI] = nominalValue ;

      const cell& c = mesh_.cells()[cellI];
      forAll(c, faceI)
      {
        const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

        if (patchID > -1 and sf.boundaryField()[patchID].size())
        {
          const scalarField& Tp(T.boundaryField()[patchID]);
          const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
          scalarField& diffP(sf.boundaryFieldRef()[patchID]);
          scalar nominalValueP;
          if (Tp[faceID] < 1673.15)
          {
            nominalValueP = 8.6e-10 * exp(-326e3/8.31446/Tp[faceID]);
          }
          else
          {
            nominalValueP = 3.7e1 * exp(-657e3/8.31446/Tp[faceID]);
          }
          diffP[faceID] = nominalValueP;
        }
      }
    }
}

void Foam::diffCoefSiC::updateAgCoef
(
  volScalarField& sf,
  const volScalarField& T,
  const labelList& addr
)
{
  forAll(addr, i)
  {
      const label cellI = addr[i];
      scalar nominalValue;
      nominalValue = 3.6e-9 * exp(-215e3/8.31446/T[cellI]);
      sf[cellI] = nominalValue ;

      const cell& c = mesh_.cells()[cellI];
      forAll(c, faceI)
      {
        const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

        if (patchID > -1 and sf.boundaryField()[patchID].size())
        {
          const scalarField& Tp(T.boundaryField()[patchID]);
          const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
          scalarField& diffP(sf.boundaryFieldRef()[patchID]);
          scalar nominalValueP;
          nominalValueP = 3.6e-9 * exp(-215e3/8.31446/Tp[faceID]);
          diffP[faceID] = nominalValueP;
        }
      }
    }
}

void Foam::diffCoefSiC::updateSrCoef
(
  volScalarField& sf,
  const volScalarField& T,
  const labelList& addr
)
{
  forAll(addr, i)
  {
      const label cellI = addr[i];
      scalar nominalValue;
      nominalValue = 1.2e-9 * exp(-205e3/8.31446/T[cellI])
                   + 1.8e6  * exp(-791e3/8.31446/T[cellI]);
      sf[cellI] = nominalValue ;

      const cell& c = mesh_.cells()[cellI];
      forAll(c, faceI)
      {
        const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

        if (patchID > -1 and sf.boundaryField()[patchID].size())
        {
          const scalarField& Tp(T.boundaryField()[patchID]);
          const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
          scalarField& diffP(sf.boundaryFieldRef()[patchID]);
          scalar nominalValueP;
          nominalValueP = 1.2e-9 * exp(-205e3/8.31446/Tp[faceID])
                        + 1.8e6  * exp(-791e3/8.31446/Tp[faceID]);
          diffP[faceID] = nominalValueP;
        }
      }
    }
}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::diffCoefSiC::updateCoef
(
    volScalarField& sf,
    const volScalarField& T,
    const labelList& addr,
    const word FpName
)
{
  //1. If the material is considered broken, then the diffCoefs for all the species
  //are set to 1e-6. //For Iaea benchmark Cases;
  if (crackFp_)
  {
    forAll(addr, i)
    {
      const label cellI = addr[i];
      scalar nominalValue;
      nominalValue = 1e-6;
      sf[cellI] = nominalValue ;

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
  //2. If the parameters are provided by users, then use them even for Cs, Ag..(input has privilege)
  forAll(FpNamesMat_, j)
  {
    if (FpNamesMat_[j] == FpName)
    {
      forAll(addr, i)
      {
        const label cellI = addr[i];
        scalar nominalValue;
        nominalValue = d1_[j] * exp(-q1_[j]/8.31446/T[cellI]) + d2_[j] * exp(-q2_[j]/8.31446/T[cellI]);
        sf[cellI] = nominalValue ;

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

  //3. If the parameters are not provided by users, then use models for Cs, Ag..
  if (FpName == "Cs")
  {
    updateCsCoef(sf, T, addr);
  }
  else if (FpName == "Kr")
  {
    updateKrCoef(sf, T, addr);
  }
  else if (FpName == "Ag")
  {
    updateAgCoef(sf, T, addr);
  }
  else if (FpName == "Sr")
  {
    updateSrCoef(sf, T, addr);
  }
  //4. If the parameters are not provided by users and we don't have a model for that
  //then error.
  else
  {
    FatalErrorIn("Foam::diffCoefSiC::updateCoef")
    << "The current model doesn't support " << FpName
    << ". Please provide its parameteres in fissionProductsTransport dictionary"
    << exit(FatalError);
  }

}





// ************************************************************************* //
