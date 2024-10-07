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

#include "poreVelocityMOXLackey.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(poreVelocityMOXLackey, 0);
    addToRunTimeSelectionTable(poreVelocityModel, poreVelocityMOXLackey, dictionary);

}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::poreVelocityMOXLackey::poreVelocityMOXLackey
(
    const fvMesh& mesh, 
    const dictionary& dict,
    const word defaultModel
)
:
    poreVelocityModel(mesh, dict, defaultModel),
    T_(nullptr),
    gradT_(nullptr),
    OM_(nullptr),
    K_(3.364)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::poreVelocityMOXLackey::~poreVelocityMOXLackey()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::poreVelocityMOXLackey::correct
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

    if( OM_ == nullptr )
    {
        OM_ = &mesh_.lookupObject<volScalarField>("oxygenMetalRatio");
    }
    
    // Reference to pore velocity field
    vectorField& poreVel = poreVelocity_.ref();

    forAll(addr, i)
    {
        const label cellI = addr[i];

        const scalar Ti = T_->internalField()[cellI];

        const vector gradTi = gradT_->internalField()[cellI];
        
        const scalar OMi = OM_->internalField()[cellI];

        const scalar logTerm = 
            - 212.275 + 65.842 * OMi + 8.9453e-2 * Ti
            - 2.55399e-2 * OMi * Ti + 2.956 * pow(OMi,2.0)
            - 5.6541e-6 * pow(Ti,2);

        const scalar p = Ti / 2000;

        // Pore velocity in cm/s
        poreVel[cellI] = 
            K_*exp(logTerm)/(p*pow(Ti,1.5))*gradTi/100.0*alphaPoreVelocity_;

        // Pore velocity in m/s
        poreVel[cellI] /= 100.0;
 
        // 
        // Correct boundary values
        // 

        // Take reference of current cell
        const cell& c = mesh_.cells()[cellI];  


        // Loop over all faces of current cell
        for(label faceI = 0 ; faceI<c.size(); ++faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

                    
            if (patchID > -1 and poreVelocity_.boundaryField()[patchID].size())
            {
                const polyPatch& pp(mesh_.boundaryMesh()[patchID]);
                const vectorField deltas = pp.faceCentres() - pp.faceCellCentres();                
                const label faceID = pp.whichFace(c[faceI]);

                // Take references boundary T
                const scalarField& Tp(T_->boundaryField()[patchID]);

                // Compute grad T at the boundary in K/cm
                const vector delta = deltas[faceID];

                const scalar magGradTp = (Tp[faceID]-Ti)/mag(delta);
                const vector gradTp = magGradTp * delta/mag(delta);

                // Take references boundary v
                vectorField& vP =
                    poreVelocity_.boundaryFieldRef()[patchID];

                // Take references boundary OM
                const scalarField& OMp(OM_->boundaryField()[patchID]);

                const scalar logTerm = 
                    - 212.275 + 65.842 * OMp[faceID] + 8.9453e-2 * Tp[faceID]
                    - 2.55399e-2 * OMp[faceID] * Tp[faceID] + 2.956 * pow(OMp[faceID],2.0)
                    - 5.6541e-6 * pow(Tp[faceID],2);

                const scalar p = Tp[faceID] / 2000;

                // cm/s
                vP[faceID] = 
                    K_*exp(logTerm)/(p*pow(Tp[faceID],1.5))*gradTp/100*alphaPoreVelocity_;

                // m/s
                vP[faceID] /= 100.0;
            }
        }
    }
}

// ************************************************************************* //
