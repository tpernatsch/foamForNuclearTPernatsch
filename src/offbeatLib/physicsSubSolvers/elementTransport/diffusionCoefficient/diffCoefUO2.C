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

#include "diffCoefUO2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diffCoefUO2, 0);
    addToRunTimeSelectionTable
    (
        diffCoefModel,
        diffCoefUO2,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffCoefUO2::diffCoefUO2
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word defaultModel
)
:
    diffCoefModel(mesh, dict, defaultModel)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffCoefUO2::~diffCoefUO2()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::diffCoefUO2::updateCsCoef
(
    volScalarField& sf,
    const volScalarField& T,
    const labelList& addr
)
{
    forAll(addr, i)
    {
        const label cellI = addr[i];
        sf[cellI] = 5.6e-8 * exp(-209e3/8.31446/T[cellI])
                  + 5.2e-4 * exp(-362e3/8.31446/T[cellI]);

        const cell& c = mesh_.cells()[cellI];
        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
            if (patchID > -1 && sf.boundaryField()[patchID].size())
            {
                const scalarField& Tp(T.boundaryField()[patchID]);
                const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                sf.boundaryFieldRef()[patchID][faceID] =
                    5.6e-8 * exp(-209e3/8.31446/Tp[faceID])
                  + 5.2e-4 * exp(-362e3/8.31446/Tp[faceID]);
            }
        }
    }
}

void Foam::diffCoefUO2::updateKrCoef
(
    volScalarField& sf,
    const volScalarField& T,
    const labelList& addr
)
{
    forAll(addr, i)
    {
        const label cellI = addr[i];
        scalar D;
        if (T[cellI] < 1773.15)
            D = 1.3e-12 * exp(-126e3/8.31446/T[cellI]);
        else
            D = 8.8e-15 * exp(-54e3/8.31446/T[cellI])
              + 6.0e-1  * exp(-480e3/8.31446/T[cellI]);
        sf[cellI] = D;

        const cell& c = mesh_.cells()[cellI];
        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
            if (patchID > -1 && sf.boundaryField()[patchID].size())
            {
                const scalarField& Tp(T.boundaryField()[patchID]);
                const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                scalar Dp;
                if (Tp[faceID] < 1773.15)
                    Dp = 1.3e-12 * exp(-126e3/8.31446/Tp[faceID]);
                else
                    Dp = 8.8e-15 * exp(-54e3/8.31446/Tp[faceID])
                       + 6.0e-1  * exp(-480e3/8.31446/Tp[faceID]);
                sf.boundaryFieldRef()[patchID][faceID] = Dp;
            }
        }
    }
}

void Foam::diffCoefUO2::updateAgCoef
(
    volScalarField& sf,
    const volScalarField& T,
    const labelList& addr
)
{
    forAll(addr, i)
    {
        const label cellI = addr[i];
        sf[cellI] = 6.7e-9 * exp(-165e3/8.31446/T[cellI]);

        const cell& c = mesh_.cells()[cellI];
        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
            if (patchID > -1 && sf.boundaryField()[patchID].size())
            {
                const scalarField& Tp(T.boundaryField()[patchID]);
                const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                sf.boundaryFieldRef()[patchID][faceID] =
                    6.7e-9 * exp(-165e3/8.31446/Tp[faceID]);
            }
        }
    }
}

void Foam::diffCoefUO2::updateSrCoef
(
    volScalarField& sf,
    const volScalarField& T,
    const labelList& addr
)
{
    forAll(addr, i)
    {
        const label cellI = addr[i];
        sf[cellI] = 2.2e-3 * exp(-488e3/8.31446/T[cellI]);

        const cell& c = mesh_.cells()[cellI];
        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
            if (patchID > -1 && sf.boundaryField()[patchID].size())
            {
                const scalarField& Tp(T.boundaryField()[patchID]);
                const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                sf.boundaryFieldRef()[patchID][faceID] =
                    2.2e-3 * exp(-488e3/8.31446/Tp[faceID]);
            }
        }
    }
}

void Foam::diffCoefUO2::updateCoef
(
    volScalarField& sf,
    const volScalarField& T,
    const labelList& addr,
    const word FpName
)
{
    // If cracked: infinite diffusivity
    if (isCracked())
    {
        forAll(addr, i)
        {
            const label cellI = addr[i];
            sf[cellI] = 1e-6;
            const cell& c = mesh_.cells()[cellI];
            forAll(c, faceI)
            {
                const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
                if (patchID > -1 && sf.boundaryField()[patchID].size())
                {
                    const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                    sf.boundaryFieldRef()[patchID][faceID] = 1e-6;
                }
            }
        }
        return;
    }

    // User-supplied per-species parameters take priority
    forAll(FpNamesMat_, j)
    {
        if (FpNamesMat_[j] == FpName)
        {
            forAll(addr, i)
            {
                const label cellI = addr[i];
                sf[cellI] = d1_[j] * exp(-q1_[j]/8.31446/T[cellI])
                          + d2_[j] * exp(-q2_[j]/8.31446/T[cellI]);
                const cell& c = mesh_.cells()[cellI];
                forAll(c, faceI)
                {
                    const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
                    if (patchID > -1 && sf.boundaryField()[patchID].size())
                    {
                        const scalarField& Tp(T.boundaryField()[patchID]);
                        const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                        sf.boundaryFieldRef()[patchID][faceID] =
                            d1_[j] * exp(-q1_[j]/8.31446/Tp[faceID])
                          + d2_[j] * exp(-q2_[j]/8.31446/Tp[faceID]);
                    }
                }
            }
            return;
        }
    }

    // Built-in correlations per species
    if (FpName == "Cs")
        updateCsCoef(sf, T, addr);
    else if (FpName == "Kr")
        updateKrCoef(sf, T, addr);
    else if (FpName == "Ag")
        updateAgCoef(sf, T, addr);
    else if (FpName == "Sr")
        updateSrCoef(sf, T, addr);
    else
        FatalErrorIn("Foam::diffCoefUO2::updateCoef")
            << "No built-in UO2 diffusion correlation for species " << FpName
            << ". Provide parameters in the diffusion dict."
            << exit(FatalError);
}

// ************************************************************************* //
