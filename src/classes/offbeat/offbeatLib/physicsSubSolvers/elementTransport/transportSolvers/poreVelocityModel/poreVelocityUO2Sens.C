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

#include "poreVelocityUO2Sens.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(poreVelocityUO2Sens, 0);
    addToRunTimeSelectionTable(poreVelocityModel, poreVelocityUO2Sens, dictionary);

}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::poreVelocityUO2Sens::poreVelocityUO2Sens
(
    const fvMesh& mesh, 
    const dictionary& dict,
    const word defaultModel
)
:
    poreVelocityModel(mesh, dict, defaultModel),
    T_(nullptr),
    gradT_(nullptr),
    R_(8.31446261815324e7)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::poreVelocityUO2Sens::~poreVelocityUO2Sens()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::poreVelocityUO2Sens::correct
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

        // Initialize correlation parameters
        scalar p0(0);
        scalar DH(0);
        if( Ti >= 1500)
        {
            p0 = exp(33.7085);
            DH = 5.98e12;
        }
        /*
        if( Ti >= 1500 and  Ti < 2000 )
        {
            p0 = exp(33.7085);
            DH = 5.98e12;
        }
        else if( Ti < 2200 )
        {
            p0 = exp(34.8107);
            DH = 6.17e12;
        }
        else if( Ti <= 2800 )
        {
            p0 = exp(35.8130);
            DH = 6.35e12;
        }
        else
        {
            p0 = exp(35.8130);
            DH = 6.35e12;

            WarningIn("void Foam::poreDiffusionSolver::updateVelocity()")
            << "Temperature is above upper bound for pore " 
            << "velocity correlation !" << endl;
        }
        */


        // Pore velocity in (cm/s)
        poreVel[cellI] =  
            alphaPoreVelocity_ * 5e-16 *
            (
                0.988 + 6.395e-6*Ti
                + 3.543e-9*pow(Ti,2.0)
                + 3e-12*pow(Ti,3.0)
            )
            * pow(Ti,-2.5) * DH * p0
            * exp(-DH/(R_*Ti)) * gradTi/100.0;

        // Convert pore velocity in m/s
        poreVel[cellI] /= 100.0;

        // 
        // Correct boundary values
        // 

        // Take reference of current cell
        const cell& c = mesh_.cells()[cellI];  

        // Loop over all faces of current cell
        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

            const polyPatch& pp(mesh_.boundaryMesh()[patchID]);

            if (patchID > -1 and poreVelocity_.boundaryField()[patchID].size())
            {
                const vectorField deltas = pp.faceCentres() - pp.faceCellCentres();
                
                const label faceID = pp.whichFace(c[faceI]);

                // Take references boundary T
                const scalarField& Tp(T_->boundaryField()[patchID]);

                // Compute grad T at the boundary
                const vector& delta = deltas[faceID];

                const scalar magGradTp = (Tp[faceID]-Ti)/mag(delta);
                const vector gradTp = magGradTp * delta/mag(delta);

                // Take references boundary v
                vectorField& vP =
                    poreVelocity_.boundaryFieldRef()[patchID];
                
                scalar p0(0);
                scalar DH(0);
                if( Tp[faceID] >= 1500 )
                {
                    p0 = exp(33.7085);
                    DH = 5.98e12;
                }
                /*
                if( Tp[faceID] >= 1500 and  Tp[faceID] < 2000 )
                {
                    p0 = exp(33.7085);
                    DH = 5.98e12;
                }
                else if( Tp[faceID] < 2200 )
                {
                    p0 = exp(34.8107);
                    DH = 6.17e12;
                }
                else if( Tp[faceID] <= 2800 )
                {
                    p0 = exp(35.8130);
                    DH = 6.35e12;
                }
                else
                {
                    p0 = exp(35.8130);
                    DH = 6.35e12;
                }
                */

                // Compute velocity in cm/s
                vP[faceID] = 
                    alphaPoreVelocity_ * 5e-16 *
                    (
                        0.988 + 6.395e-6*Tp[faceID]
                        + 3.543e-9*pow(Tp[faceID],2.0)
                        + 3e-12*pow(Tp[faceID],3.0)
                    )
                    * pow(Tp[faceID],-2.5) * DH * p0
                    * exp(-DH/(R_*Tp[faceID]))
                    * gradTp/100;

                // Convert velocity to m/s
                vP[faceID] /= 100.0;

            }
        }
    }
}

// ************************************************************************* //
