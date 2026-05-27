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

#include "poreVelocityMOXClementFinnis.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(poreVelocityMOXClementFinnis, 0);
    addToRunTimeSelectionTable(poreVelocityModel, poreVelocityMOXClementFinnis, dictionary);

}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::poreVelocityMOXClementFinnis::poreVelocityMOXClementFinnis
(
    const fvMesh& mesh, 
    const dictionary& dict,
    const word defaultModel
)
:
    poreVelocityModel(mesh, dict, defaultModel),
    T_(nullptr),
    gradT_(nullptr),
    A_(4220),
    alpha_(66490),
    l_(80e-6)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::poreVelocityMOXClementFinnis::~poreVelocityMOXClementFinnis()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::poreVelocityMOXClementFinnis::correct
(
    const labelList& addr
)
{
    if( T_ == nullptr )
    {
        T_ = &mesh_.lookupObject<volScalarField>("T");
    }

    if( gradT_ == nullptr )
    {
        gradT_ = &mesh_.lookupObject<volVectorField>("gradT");
    }

    // Reference to pore velocity field
    vectorField& poreVel = poreVelocity_.ref();
    
    forAll(addr, i)
    {
        const label cellI = addr[i];

        const scalar Ti = T_->internalField()[cellI];

        const vector gradTi = gradT_->internalField()[cellI];

        // Pore velocity in m/s
        poreVel = A_*exp(-alpha_/Ti)*l_*gradTi*alphaPoreVelocity_;
 
        // 
        // Correct boundary values
        // 

        // Take reference of current cell
        const cell& c = mesh_.cells()[cellI];  

        // Loop over all faces of current cell
        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

            if (patchID > -1 and poreVelocity_.boundaryField()[patchID].size())
            {
                const polyPatch& pp(mesh_.boundaryMesh()[patchID]);
                
                const vectorField deltas = pp.faceCentres() - pp.faceCellCentres();
                
                const label faceID = pp.whichFace(c[faceI]);

                // Take references boundary T
                const scalarField& Tp(T_->boundaryField()[patchID]);

                // Compute grad T at the boundary in K/m
                const vector delta = deltas[faceID];

                const scalar magGradTp = (Tp[faceID]-Ti)/mag(delta);
                const vector gradTp = magGradTp * delta/mag(delta);

                // Take references boundary v
                vectorField& vP =
                    poreVelocity_.boundaryFieldRef()[patchID];

                vP[faceID] = A_*exp(-alpha_/Tp[faceID])*l_*gradTp*alphaPoreVelocity_;
            }
        }
    }
}

// ************************************************************************* //
