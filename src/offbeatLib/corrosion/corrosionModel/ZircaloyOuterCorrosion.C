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

#include "ZircaloyOuterCorrosion.H"
#include "addToRunTimeSelectionTable.H"
#include "resistiveLayerFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ZircaloyOuterCorrosion, 0);

    addToRunTimeSelectionTable
    (
        corrosionModel,
        ZircaloyOuterCorrosion,
        dictionary
    );
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ZircaloyOuterCorrosion::ZircaloyOuterCorrosion
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    corrosionModel(mesh, dict),
    kineticsModel_
    (
        oxidationKineticsModel::New
        (
            word(dict.lookup("oxidationKineticsModel")),
            dict
        )
    )
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ZircaloyOuterCorrosion::~ZircaloyOuterCorrosion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::ZircaloyOuterCorrosion::correct
(
    surfaceScalarField& oxideThickness,
    surfaceScalarField& DOxideThickness,
    const label& patchID
)
{
    // Legend for the following lines:
    //  - Tp cell-center temperature for cells attached to patch
    //  - Tb face-center temperature on outer surface of oxide
    //  - Ti - temperature at oxide/metal interface

    // Current and oldTime oxide layer thickness
    const fvPatch& patch(mesh_.boundary()[patchID]);
    scalarField& S(oxideThickness.boundaryFieldRef()[patchID]);
    const scalarField& Sold(oxideThickness.oldTime().boundaryField()[patchID]);

    // Oxide layer increase
    scalarField& DS(DOxideThickness.boundaryFieldRef()[patchID]);

    // Delta time
    scalar deltaT = mesh_.time().deltaTValue();

    // Inverse cell to patch face distances
    // Inverse cell to patch face distances
    const scalarField d = patch.deltaCoeffs();

    // Temperature on the patch (= metal/oxide interface temperature)
    const fvPatchScalarField& Tb
        = patch.lookupPatchField<volScalarField, scalar>("T");

    // Find temperatures of cells adjacents to patch
    const scalarField Tp = Tb.patchInternalField();

    // Fast flux on the patch
    const scalarField& phi
        = patch.lookupPatchField<volScalarField, scalar>("fastFlux");

    // Thermal conductivity
    fvPatchScalarField& kb = const_cast<fvPatchScalarField&>
    (
        patch.lookupPatchField<volScalarField, scalar>("k")
    );

    // Pilling-Bedworth ratio for Zr
    const scalar Rpb = 1.56;

    const scalarField kp = kb.patchInternalField();
    scalarField Ti(Tb.size());
    scalarField kox(Tb.size());
    scalarField alpha_ox(Tb.size());

    // Calculate the boundary effective thermal conductivity and
    // interface temperature
    forAll(patch, i)
    {
        // Oxide thermal properties
        // Dollins (1983) correlation
        kox[i] = 1.9599-(2.41e-4 - (6.43e-7 - 1.94e-10*Tb[i])*Tb[i])*Tb[i];
        
        // MATPRO correlation
        //kox[i] = 0.835 + 0.000181*Tb[i];
        const scalar dox = max(0.5*(Sold[i] + S[i]), SMALL);

        // Eroded metal thickness
        const scalar dmetal = dox/Rpb;

        // Effective HTC of oxide layer minus the eroded metal
        alpha_ox[i] = 1/(dox/kox[i] - dmetal/kp[i]);

        // Effective HTC
        const scalar alpha_m = kp[i]*d[i];
        const scalar beta = alpha_ox[i] / (alpha_ox[i] + alpha_m);

        //- Correct the boundary conductivity
        kb[i] = beta*kp[i];

        // Interface temperature
        Ti[i] = beta*Tb[i] + (1-beta)*Tp[i];
    }
    
    if (debug)
    {
        const scalarField& Ap = patch.magSf();
        scalar Atot = gSum(Ap);
        
        Info<< "Corrosion for patch " << patch.name() << nl
            << tab << "Internal temperature: " << gSum(Tp*Ap)/Atot << nl
            << tab << "Interface temperature: " << gSum(Ti*Ap)/Atot << nl
            << tab << "Boundary temperature: " << gSum(Tb*Ap)/Atot << nl
            << tab << "Internal conductivity: " << gSum(kp*Ap)/Atot << nl
            << tab << "Effective conductivity: " << gSum(kb*Ap)/Atot << nl
            << tab << "Oxide thermal conductivity: " << gSum(kox*Ap)/Atot << nl
            << endl;
    }

    if (isA<resistiveLayerFvPatchScalarField>(kb))
    {
        resistiveLayerFvPatchScalarField& pf = refCast<resistiveLayerFvPatchScalarField>(kb);
        pf.alpha() = alpha_ox;
    }
    else
    {
        FatalErrorInFunction()
            << "Patch " << patch.name() << " on field "
            << kb.internalField().name() << " must be of type "
            << resistiveLayerFvPatchScalarField::typeName
            << " to model the oxide interface temperature"
            << " correctly. Current type is " << kb.type()
            << "." << endl << abort(FatalError);
    }

    kineticsModel_->correctOxideThickness(S, Sold, Ti, Tb, phi, deltaT);

    DS = S - Sold;

    updateMesh_ = false;
};


void Foam::ZircaloyOuterCorrosion::updateDMetalThickness
(
    const surfaceScalarField& DOxideThickness,
    surfaceScalarField& DMetalThickness,
    const label& patchID
)
{
    // Pilling-Bedworth ratio for Zr
    const scalar Rpb = 1.56;

    // Oxide layer increase
    const scalarField& DS_oxide(DOxideThickness.boundaryField()[patchID]);
    scalarField& DS_metal =(DMetalThickness.boundaryFieldRef()[patchID]);

    DS_metal = DS_oxide/Rpb;
};

// ************************************************************************* //
