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

#include "diffCoefPiecewiseArrhenius.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(diffCoefPiecewiseArrhenius, 0);
    addToRunTimeSelectionTable(diffCoefModel, diffCoefPiecewiseArrhenius, dictionary);

    static const scalar R = 8.31446;   // J/(mol·K)
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffCoefPiecewiseArrhenius::diffCoefPiecewiseArrhenius
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
    q2_(readScalar(dict.lookup("q2"))),
    Tswitch_(readScalar(dict.lookup("Tswitch")))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

namespace
{
    inline Foam::scalar piecewise
    (
        Foam::scalar T,
        Foam::scalar d1, Foam::scalar q1,
        Foam::scalar d2, Foam::scalar q2,
        Foam::scalar Tswitch, Foam::scalar R
    )
    {
        return T < Tswitch
            ? d1*exp(-q1/R/T)
            : d2*exp(-q2/R/T);
    }
}

void Foam::diffCoefPiecewiseArrhenius::updateCoef
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
        sf[cellI] = piecewise(T[cellI], d1_, q1_, d2_, q2_, Tswitch_, R);

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
                    piecewise(Tp[faceID], d1_, q1_, d2_, q2_, Tswitch_, R);
            }
        }
    }
}

// ************************************************************************* //
