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
#include "temperatureCoupledFvPatchScalarField.H"
#include "Time.H"
#include "gapGasModel.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Private members  * * * * * * * * * * * * *//

Foam::tmp<Foam::scalarField>
Foam::temperatureCoupledFvPatchScalarField::weights() const
{

    if (coupled_)
    {
        const fvPatch& patch = this->patch();
        const fvPatch& nbrPatch = regionCoupledPatch_.neighbFvPatch();

        const scalarField deltas
        (
            patch.nf() & patch.delta()
        );

        const scalarField alphaP(kappa()/deltas);

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

        const scalarField alphaN(nbrKappa/max(nbrDeltas, SMALL));       

        tmp<scalarField> tw(new scalarField(deltas.size()));
        scalarField& w = tw.ref();

        forAll(alphaP, facei)
        {
            scalar di = alphaP[facei];
            scalar dni = alphaN[facei];

            w[facei] = di/(di + dni);
        }

        return tw;
    }
    else
    {
        return patch().weights();
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperatureCoupledFvPatchScalarField::
temperatureCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCoupledFvPatchScalarField(p, iF),
    kappaName_("unknown")
{
    regionCoupledFvPatchScalarField::coupled_ = true;
}


Foam::temperatureCoupledFvPatchScalarField::
temperatureCoupledFvPatchScalarField
(
    const temperatureCoupledFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    regionCoupledFvPatchScalarField(ptf, p, iF, mapper),
    kappaName_(ptf.kappaName_)
{
    regionCoupledFvPatchScalarField::coupled_ = true;
}


Foam::temperatureCoupledFvPatchScalarField::
temperatureCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    regionCoupledFvPatchScalarField(p, iF, dict),
    kappaName_(dict.lookupOrDefault<word>("kappa", "k"))
{
    regionCoupledFvPatchScalarField::coupled_ = true;
}


Foam::temperatureCoupledFvPatchScalarField::
temperatureCoupledFvPatchScalarField
(
    const temperatureCoupledFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCoupledFvPatchScalarField(ptf, iF),
    kappaName_(ptf.kappaName_)
{
    regionCoupledFvPatchScalarField::coupled_ = true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::temperatureCoupledFvPatchScalarField::write(Ostream& os) const
{
    regionCoupledFvPatchScalarField::write(os);

#ifdef OPENFOAMFOUNDATION        
    writeEntryIfDifferent<word>(os, "kappa", "k", kappaName_);
#elif OPENFOAMESI
    os.writeEntryIfDifferent<word>("kappa", "k", kappaName_);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        temperatureCoupledFvPatchScalarField
    );
};


// ************************************************************************* //
