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

#include "TSSpDXD4Liner.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TSSpDXD4Liner, 0);
    addToRunTimeSelectionTable
    (
        TSSpModel,
        TSSpDXD4Liner, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TSSpDXD4Liner::TSSpDXD4Liner
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& lawDict
)
:
    TSSpModel(mesh, mat, lawDict)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TSSpDXD4Liner::~TSSpDXD4Liner()
{}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::TSSpDXD4Liner::correctTSSp
(
    const labelList& addr
)
{
    const volScalarField& T_(mesh_.lookupObject<volScalarField>("T"));
    const scalarField& Ti = T_.internalField();

    //TODO: Should be able to read the TSSp0 and Qp for the selected TSSP model for substrate
    //Using Kammenzind's expression for TSSp substrate
    const scalar TSSp0 = 31000;
    const scalar Qp = 25239;

    //TODO: Should be read in the solverDict if TSSp for liner is chosen as TSSpDDXD4Liner
    //chemical potential difference from experimental results of Wong et al. (Uncertainty of ~42%)
    const scalar mu_H0 = -656; // [J/mol]

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        // TSSp_liner = TSSp_substrate * exp(mu_H0/RT) {From Gong et al.}
        TSSp_[cellI] = (TSSp0 * exp(-Qp / (8.31446 * Ti[cellI]))) * exp(mu_H0 / (8.31446 * Ti[cellI]));

        const cell& c = mesh_.cells()[cellI];  
        forAll(c, faceI)
        {
            const label patchID = 
            mesh_.boundaryMesh().whichPatch(c[faceI]);

            if (patchID > -1 and TSSp_.boundaryField()[patchID].size())
            {
                const label faceID = 
                mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                scalarField& TSSpP(TSSp_.boundaryFieldRef()[patchID]);
                    
                const scalarField& Tp(T_.boundaryField()[patchID]);
                TSSpP[faceID] = (TSSp0 * exp(-Qp / (8.31446 * Tp[faceID]))) * exp(mu_H0 / (8.31446 * Tp[faceID]));
            }
        }
    }
    TSSp_.correctBoundaryConditions();
}
// ************************************************************************* //
