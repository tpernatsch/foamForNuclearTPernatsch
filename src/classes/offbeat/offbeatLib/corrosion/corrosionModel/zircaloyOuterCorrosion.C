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

#include "zircaloyOuterCorrosion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zircaloyOuterCorrosion, 0);

    addToRunTimeSelectionTable
    (
        corrosionModel, 
        zircaloyOuterCorrosion, 
        dictionary
    );
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zircaloyOuterCorrosion::zircaloyOuterCorrosion
(
    const fvMesh& mesh
)
:
    corrosionModel(mesh)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zircaloyOuterCorrosion::~zircaloyOuterCorrosion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::zircaloyOuterCorrosion::correct
(
    surfaceScalarField& oxideThickness,
    surfaceScalarField& DOxideThickness,
    const label& patchID
)
{
    // Legend for the following lines:
    //  - Ti cell-center temperature for cells attached to patch
    //  - Tf face-center temperature, interface between metal and oxide
    //  - To outer face-center temperature, oxide surface

    const fvPatch& patch(mesh_.boundary()[patchID]);

    // Current and oldTime oxide layer thickness
    scalarField& S(oxideThickness.boundaryFieldRef()[patchID]);
    const scalarField& Sold(oxideThickness.oldTime().boundaryField()[patchID]);
    const scalarField& Sprev(oxideThickness.prevIter().boundaryField()[patchID]);

    // Oxide layer increase 
    scalarField& DS(DOxideThickness.boundaryFieldRef()[patchID]);

    // Delta time
    scalar deltaT = mesh_.time().deltaTValue();

    // Inverse cell to patch face distances
    const scalarField d = patch.deltaCoeffs();

    // Temperature on the patch (= metal/oxide interface temperature)
    const scalarField& Tf
    = patch.lookupPatchField<volScalarField, scalar>("T");

    // Find temperatures of cells adjacents to patch 
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    const labelList faceCells = patch.faceCells();
    const scalarField Ti(T.internalField(), faceCells);

    // Fast flux on the patch
    const scalarField& phi 
    = patch.lookupPatchField<volScalarField, scalar>("fastFlux");

    // Loop over all the faces of the boundary patch
    forAll(patch, i)
    {    
        // Normal operating condition model (EPRI KWU CE) 
        if ( Tf[i] < 673 )
        {
            // Correlation parameters
            const scalar C1(6.3e-9);
            const scalar C2((8.04e7+2.59*1e8*pow(7.46e-15*phi[i],0.25))*1e-6);
            const scalar Q1R(16266);
            const scalar Q2R(13775);

            // Transition thickness
            scalar sTransition = 2e-6 ;

            // Compute oxide thickness
            if( S[i] <= sTransition )
            {
                S[i] = pow
                (
                    C1*exp(-Q1R/Tf[i])*deltaT/3600.0/24.0
                    + pow(Sold[i], 3.0)
                    , 1.0/3.0
                );
            }
            // If S is crossing sTransition in this timeStep
            else if( Sold[i] < sTransition )
            {
                // Calculate S continuing with pre-transition law
                scalar S_preTransition = pow
                (
                    C1*exp(-Q1R/Tf[i])*deltaT/3600.0/24.0
                    + pow(Sold[i], 3.0)
                    , 1.0/3.0
                );

                // Calculate the deltaT up to sTransition
                scalar deltaTpreTransition = deltaT*
                (sTransition - Sold[i])/(S_preTransition - Sold[i]);

                // Calculate the deltaT post sTransition
                scalar deltaTpostTransition = deltaT - deltaTpreTransition;

                S[i] = sTransition 
                + C2*exp(-Q2R/Tf[i])*deltaTpostTransition/3600.0/24.0;
            }
            else
            {
                S[i] = Sold[i] + C2*exp(-Q2R/Tf[i])*deltaT/3600.0/24.0;
            }
        }
        // High temperature model:
        else
        {
            // Correlation parameters
            scalar As;
            scalar QsR;

            // Leistikov correlation
            if ( Tf[i] < 1800 )
            {
                As  = 7.82e-6;
                QsR = 20214;
            }
            // Interpolation (Procedure from G. Schanz - 2003)
            else if (Tf[i] < 1900)
            {
                QsR = 1900*1800/(1900-1800)
                    * log
                      (
                        2.98e-3*Foam::exp(-28420.0/1900.0)
                        /
                        7.82e-6*Foam::exp(-20214.0/1800.0)
                      );

                As  = 7.82e-6*exp(QsR/1900.0);
            }
            // Prater-Courtright correlation
            else
            {
                As  = 2.98e-3;
                QsR = 28420;
            }

            S[i] = sqrt( As*exp(-QsR/Tf[i])*deltaT + pow(Sold[i], 2.0) );
        }

        DS[i] = S[i] - Sold[i];
    }

    updateMesh_ = false;

    // Check absolute change and decide whether to update the mesh or not
    scalar oxideChange(gMax(mag(S - Sprev)));
    if(oxideChange > 1e-10)
    {
        updateMesh_ = true;
    }
};


void Foam::zircaloyOuterCorrosion::updateDMetalThickness
(
    const surfaceScalarField& DOxideThickness,
    surfaceScalarField& DMetalThickness,
    const label& patchID
)
{
    // Pilling-Bedworth ratio for Zr
    scalar Rpb(1.56);

    // Oxide layer increase 
    const scalarField& DS_oxide(DOxideThickness.boundaryField()[patchID]);

    // Change in metal layer
    scalarField& DS_metal =(DMetalThickness.boundaryFieldRef()[patchID]);

    DS_metal = DS_oxide/Rpb;
};

// ************************************************************************* //
