/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "oxidePickupFractionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

oxidePickupFractionFvPatchScalarField::oxidePickupFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    hydrogenTransportFvPatchScalarField(p, iF),
    MH_(1.00784),
    MZr_(91.224),
    volFactorFlag_(false),
    rInner_(4.565),
    rOuter_(5.315),
    pickupFrac_(0.15)
{}

oxidePickupFractionFvPatchScalarField::oxidePickupFractionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    hydrogenTransportFvPatchScalarField(p, iF, dict),
    MH_(1.00784),
    MZr_(91.224),
    volFactorFlag_(dict.lookupOrDefault<bool>("volFactor", false)),
    rInner_(4.5),
    rOuter_(5.32),
    pickupFrac_(dict.lookupOrDefault<scalar>("pickupFraction", 0.15))
{
    if (volFactorFlag_)
    {
    #ifdef OPENFOAMFOUNDATION
        rInner_ = dict.lookup<scalar>("rInner");
        rOuter_ = dict.lookup<scalar>("rOuter");
    #elif OPENFOAMESI
        rInner_ = dict.get<scalar>("rInner");
        rOuter_ = dict.get<scalar>("rOuter");
    #endif
    }
}

oxidePickupFractionFvPatchScalarField::oxidePickupFractionFvPatchScalarField
(
    const oxidePickupFractionFvPatchScalarField& tdpvf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    hydrogenTransportFvPatchScalarField(tdpvf, p, iF, mapper),
    MH_(1.00784),
    MZr_(91.224),
    volFactorFlag_(false),
    rInner_(4.565),
    rOuter_(5.315),
    pickupFrac_(0.15)
{}


oxidePickupFractionFvPatchScalarField::oxidePickupFractionFvPatchScalarField
(
    const oxidePickupFractionFvPatchScalarField& tdpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    hydrogenTransportFvPatchScalarField(tdpvf, iF),
    MH_(tdpvf.MH_),
    MZr_(tdpvf.MZr_),
    volFactorFlag_(tdpvf.volFactorFlag_),
    rInner_(tdpvf.rInner_),
    rOuter_(tdpvf.rOuter_),
    pickupFrac_(tdpvf.pickupFrac_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void oxidePickupFractionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    hydrogenTransportFvPatchScalarField::autoMap(m);
}


void oxidePickupFractionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    hydrogenTransportFvPatchScalarField::rmap(ptf, addr);
}


void oxidePickupFractionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //- Entering hydrogen flux from oxydation rate of the clad outer wall
    // Reference to dOxide patch field
    const scalarField& dOx =
        this->patch().lookupPatchField<surfaceScalarField, scalar>("DOxideThickness");
    // Reference to Oxide patch field
    const scalarField& oxTh =
        this->patch().lookupPatchField<surfaceScalarField, scalar>("oxideThickness");

    //- Oxide thickness growthRate
    const scalar dt = this->db().time().deltaTValue();
    scalarField growthRate = dOx/dt;

    // Ingress rate from Courty, Motta, et al.
    // Pilling-Bedworth ratio for Zr = 1.56
    ingressRate_ = 1e6 * 4 * pickupFrac_ / 1.56 * MH_ / MZr_ * growthRate;

    if(volFactorFlag_)
    {
        //- Oxide thickness average
        scalarField oxAvg = (oxTh + (oxTh - dOx)) / 2;
        scalar cladTh = rOuter_ - rInner_;

        scalarField volFactor = (2*rOuter_*oxAvg - pow(oxAvg, 2)) /
            (2*rOuter_ * (cladTh-oxAvg) - pow(cladTh, 2) + pow(oxAvg, 2));

        ingressRate_ *= volFactor;
    }

    if (debug)
    {
        Info<< "Hydrogen uptake summary: " << nl
            << tab << "thickness: " << oxTh << nl
            << tab << "dOx: " << dOx << nl
            << tab << "growth rate: " << growthRate << nl
            << tab << "Ingress rate: " << ingressRate_ << nl
            << endl;
    }

    hydrogenTransportFvPatchScalarField::updateCoeffs();
}

void oxidePickupFractionFvPatchScalarField::write(Ostream& os) const
{
    hydrogenTransportFvPatchScalarField::write(os);
#ifdef OPENFOAMFOUNDATION
    if (volFactorFlag_)
    {
        writeEntryIfDifferent<bool>(os, "volFactor", false, volFactorFlag_);
        writeEntryIfDifferent<scalar>(os, "rInner", 4.565, rInner_);
        writeEntryIfDifferent<scalar>(os, "rOuter", 5.315, rOuter_);
    }

    writeEntryIfDifferent<scalar>(os, "pickupFraction", 0.15, pickupFrac_);

#elif OPENFOAMESI
    if (volFactorFlag_)
    {
        os.writeEntryIfDifferent<bool>("volFactor", false, volFactorFlag_);
        os.writeEntryIfDifferent<scalar>("rInner", 4.565, rInner_);
        os.writeEntryIfDifferent<scalar>("rOuter", 5.315, rOuter_);
    }

    os.writeEntryIfDifferent<scalar>("pickupFraction", 0.15, pickupFrac_);
#endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    oxidePickupFractionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
