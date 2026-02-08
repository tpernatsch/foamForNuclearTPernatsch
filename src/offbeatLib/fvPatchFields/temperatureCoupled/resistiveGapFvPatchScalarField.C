/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"
#include "resistiveGapFvPatchScalarField.H"
#include "Time.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Private members  * * * * * * * * * * * * *//


Foam::tmp<Foam::scalarField>
Foam::resistiveGapFvPatchScalarField::weights() const
{
    if (coupled_)
    {
        const fvPatch& patch = this->patch();
        const fvPatch& nbrPatch = regionCoupledPatch_.neighbFvPatch();

        const scalarField deltas
        (
            patch.nf() & patch.delta()
        );

        const scalarField nbrKappa
        (
            regionCoupledPatch_.regionCoupledPatch().interpolate
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>(kappaName_)
            )
        );

        const scalarField nbrDeltas
        (
            regionCoupledPatch_.regionCoupledPatch().interpolate
            (
                nbrPatch.nf() & nbrPatch.delta()
            )
        );

        const scalarField alphaP(kappa()/deltas);

        const scalarField alphaN(nbrKappa/max(nbrDeltas, SMALL));
        
        scalarField alphaEff = 1/(1/alphaP + 1/max(alphaN, SMALL) + 1/max(alpha_gap_, SMALL));
        
        return 1 - alphaEff / alphaP;
    }
    else
    {
        return patch().weights();
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::resistiveGapFvPatchScalarField::
resistiveGapFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    temperatureCoupledFvPatchScalarField(p, iF),
    alpha_gap_(this->size(), 0)
{}


Foam::resistiveGapFvPatchScalarField::
resistiveGapFvPatchScalarField
(
    const resistiveGapFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    temperatureCoupledFvPatchScalarField(ptf, p, iF, mapper),
    alpha_gap_(mapper(ptf.alpha_gap_))
{}


Foam::resistiveGapFvPatchScalarField::
resistiveGapFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    temperatureCoupledFvPatchScalarField(p, iF, dict),
    alpha_gap_(scalarField("alpha", dict, p.size()))
{}


Foam::resistiveGapFvPatchScalarField::
resistiveGapFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const bool alphaRequired
)
:
    temperatureCoupledFvPatchScalarField(p, iF, dict),
    alpha_gap_(p.size(), 0)
{
    if( alphaRequired )
    {
        alpha_gap_ = scalarField("alpha", dict, p.size());
    }
}

Foam::resistiveGapFvPatchScalarField::
resistiveGapFvPatchScalarField
(
    const resistiveGapFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    temperatureCoupledFvPatchScalarField(ptf, iF),
    alpha_gap_(ptf.alpha_gap_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::resistiveGapFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    temperatureCoupledFvPatchScalarField::autoMap(m);
#ifdef OPENFOAMFOUNDATION    
    m(alpha_gap_, alpha_gap_);
#elif OPENFOAMESI
    alpha_gap_.autoMap(m);
#endif    
}


void Foam::resistiveGapFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    temperatureCoupledFvPatchScalarField::rmap(ptf, addr);

    const resistiveGapFvPatchScalarField& dmptf =
        refCast<const resistiveGapFvPatchScalarField>(ptf);

    alpha_gap_.rmap(dmptf.alpha_gap_, addr);
}


void Foam::resistiveGapFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    if (!regionCoupledPatch_.owner())
    {
        // Update the heat transfer coefficient from the owner patchField
        const fvPatch& nbrPatch = regionCoupledPatch_.neighbFvPatch();
        
        alpha_gap_ = regionCoupledPatch_.regionCoupledPatch().interpolate
        (
            refCast<const resistiveGapFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>
                (
                    internalField().name()
                )
            ).alpha()
        );
    }

    const_cast<fvPatchScalarField&>
    (
        patch().lookupPatchField<volScalarField, scalar>("hGap")
    ) = alpha_gap_;
    
    // The if statement should be performed for the same patch (owner or nbr)
    // where the hGap is directly updated (not mapped from the other side)
    if (regionCoupledPatch_.owner())
    {
        if( gMax(alpha_gap_) <= 0 )
        {
            FatalErrorInFunction
            << patch().name() << " has no positive gap conductance." 
            << exit(FatalError);
        }
    }

    temperatureCoupledFvPatchScalarField::updateCoeffs();
}


void Foam::resistiveGapFvPatchScalarField::write(Ostream& os) const
{
    temperatureCoupledFvPatchScalarField::write(os);
    
#ifdef OPENFOAMFOUNDATION        
    writeEntry(os, "alpha", alpha_gap_);
#elif OPENFOAMESI
    alpha_gap_.writeEntry("alpha", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        resistiveGapFvPatchScalarField
    );
};


// ************************************************************************* //
