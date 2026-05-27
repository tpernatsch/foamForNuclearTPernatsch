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

#include "diffCoefDoubleArrhenius.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(diffCoefDoubleArrhenius, 0);
    addToRunTimeSelectionTable(diffCoefModel, diffCoefDoubleArrhenius, dictionary);

    static const scalar R = 8.31446;   // J/(mol·K)
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffCoefDoubleArrhenius::diffCoefDoubleArrhenius
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word defaultModel
)
:
    diffCoefModel(mesh, dict, defaultModel),
    d1_(readScalar(dict.lookup("d1"))),
    q1_(readScalar(dict.lookup("q1"))),
    d2_(readScalar(dict.lookup("d2"))),
    q2_(readScalar(dict.lookup("q2")))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::diffCoefDoubleArrhenius::updateCoef
(
    volScalarField& sf,
    const volScalarField& T,
    const labelList& addr,
    const word /*FpName*/
)
{
    forAll(addr, i)
    {
        const label cellI = addr[i];
        sf[cellI] = d1_*exp(-q1_/R/T[cellI]) + d2_*exp(-q2_/R/T[cellI]);

        const cell& c = mesh_.cells()[cellI];
        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
            if (patchID > -1 && sf.boundaryField()[patchID].size())
            {
                const label faceID =
                    mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                const scalarField& Tp = T.boundaryField()[patchID];
                sf.boundaryFieldRef()[patchID][faceID] =
                    d1_*exp(-q1_/R/Tp[faceID]) + d2_*exp(-q2_/R/Tp[faceID]);
            }
        }
    }
}

// ************************************************************************* //
