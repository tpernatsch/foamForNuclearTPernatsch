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
#include "trisoGapFvPatchScalarField.H"
#include "regionCoupledOFFBEATFvPatch.H"
#include "gapGasModel.H"
#include "physicoChemicalConstants.H"
#include "gapContactFvPatchVectorField.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Private members  * * * * * * * * * * * * *//


Foam::tmp<Foam::scalarField>
Foam::trisoGapFvPatchScalarField::gapWidth() const
{

    const regionCoupledOFFBEATFvPatch& patch
        = refCast<const regionCoupledOFFBEATFvPatch>(regionCoupledPatch_);

    const regionCoupledOFFBEATFvPatch& nbrPatch
        = refCast<const regionCoupledOFFBEATFvPatch>(patch.nbrPatch());

    const fvPatchVectorField& totalDispPatch =
    (
        patch.boundaryMesh().mesh().foundObject<fvMesh>("referenceMesh")
    ) ?
    patch.lookupPatchField<volVectorField, vector>("DD")
    :
    patch.lookupPatchField<volVectorField, vector>("D");

    const fvPatchVectorField& totalDispNbrPatch =
    (
        patch.boundaryMesh().mesh().foundObject<fvMesh>("referenceMesh")
    ) ?
    nbrPatch.lookupPatchField<volVectorField, vector>("DD")
    :
    nbrPatch.lookupPatchField<volVectorField, vector>("D");

    vectorField nf = patch.Sf() / patch.magSf();
    vectorField nbrNf = -regionCoupledPatch_.regionCoupledPatch().interpolate
                      (
                             nbrPatch.Sf()/nbrPatch.magSf()
                      );

    vectorField Cf = patch.Cf()
                   + totalDispPatch;

    vectorField nbrCf = regionCoupledPatch_.regionCoupledPatch().interpolate
                      (
                             nbrPatch.Cf()
                           + totalDispNbrPatch
                      );

    //- The radii of current patch and its neighbour patch.
    r1_ =  mag(Cf & nf);
    r2_ = mag(nbrCf & nf);


    return max((nbrCf - Cf) & nf, 0.0);
}


Foam::tmp<Foam::scalarField>
Foam::trisoGapFvPatchScalarField::hGap() const
{
    const gapGasModel& gapGas = patch().boundaryMesh().mesh()
            .lookupObject<gapGasModel>(gapGasModelName_);

    const regionCoupledOFFBEATFvPatch& patch
        = refCast<const regionCoupledOFFBEATFvPatch>(regionCoupledPatch_);

    const regionCoupledOFFBEATFvPatch& nbrPatch
        = refCast<const regionCoupledOFFBEATFvPatch>(patch.nbrPatch());

    const trisoGapFvPatchScalarField& nbr = this->neighbour();

    const scalarField& pT = *this;
    scalarField nbrT = patch.regionCoupledPatch().interpolate(nbr);
    scalarField nbrR = patch.regionCoupledPatch().interpolate(nbr.R());

    scalarField nbrjd = patch.regionCoupledPatch().interpolate(nbr.jd());

    const scalarField& k
        = patch.lookupPatchField<volScalarField, scalar>(kappaName_);

    scalarField nbrK = patch.regionCoupledPatch().interpolate
    (
        nbrPatch.lookupPatchField<volScalarField, scalar>(kappaName_)
    );

    const scalarField& eps
        = patch.lookupPatchField<volScalarField, scalar>(emissivityName_);

    scalarField nbrEps = patch.regionCoupledPatch().interpolate
    (
        nbrPatch.lookupPatchField<volScalarField, scalar>(emissivityName_)
    );

    const scalarField& interfaceP(patch.lookupPatchField<volScalarField, scalar>("interfaceP"));

    scalarField hGas(this->size(), 0);
    scalarField hContact(this->size(), 0);
    scalarField hRad(this->size(), 0);

    forAll(pT,  i)
    {
        // Gas HTC
        //- The gas conduction for a sphere have a form like:
        // hGas = k / denominator
        // denominator = r²(1/rIn - 1/rOut) + Cr(R1 + R2) + g1 + rIn²/rOut²g2;
        // Cr​ is a roughness coefficient with R1, R2 as the surface roughness, we use
        // FRAPCON model to calculate this term - Cr(R1 + R2);
        // rIn, rOut is the inner and outer radius of the gap, r is the radius of
        // current patch.
        // g1 and g2 are the jump distance for the 2 patches. The current model only
        // allow us to have the value of (g1 + g2), so instead of g1 + rIn²/rOut²g2
        // we use the FRAPCON model(g1 + g2), approximatly
        {
            // Average interface temperature
            scalar T = 0.5*(pT[i] + nbrT[i]);
            scalar d_jump = 0.0;

            // Gap gas conductivity
            const scalar k = gapGas.kappa(T);

            // Interface pressure assumed zero (i.e. open gap)
            scalar pI = interfaceP[i]/1e4/9.8;

            // Gap width
            scalar rIn = min(r1_[i], r2_[i]);
            scalar rOut = max(r1_[i], r2_[i]);
            scalar gap_sphere = max(r1_[i]*r1_[i]*(1.0/rIn - 1.0/rOut), 0.0);

            // roughness distance
            scalar d_eff = 0;

            d_eff = exp(-1.25e-3*pI)*(R_[i] + nbrR[i]);

            // If the jump distance is provided, use the provided value.
            if (jd_[i] < 0.0 && nbrjd[i] < 0.0 )
            {
              d_jump = 0.0137*k*sqrt(T) / max((gapGas.p()*gapGas.a(T)),SMALL);
            }
            else if (jd_[i] >=0.0 && nbrjd[i] >= 0.0)
            {
              d_jump = jd_[i] + nbrjd[i];
            }
            else
            {
              FatalErrorIn
              (
                  "Foam::trisoGapFvPatchScalarField::hGap() "
              ) << "If the jumpDistance is provided, it must be provided for both patches"
                << exit(FatalError);
            }

            scalar d_gap_total = gap_sphere + d_eff;

            hGas[i] = k / max((d_gap_total + d_jump), SMALL);
            // hGas[i] = k / (d_gap + 1.8*d_jump + d_eff );
            hGas[i] = max(hGas[i], 0.0);
        }
        // Radiative HTC
        {
            const scalar T1 = pT[i];
            const scalar T2 = nbrT[i];
            const scalar eps1 = max(eps[i], SMALL);
            const scalar eps2 = max(nbrEps[i], SMALL);

            hRad[i] = constant::physicoChemical::sigma.value()*
                      (T1 + T2)*(T1*T1 + T2*T2)/(1/eps1 + 1/eps2 - 1);

        }
        // Contact HTC
        {
            const scalar meyerHardness = 680e6;

            scalar Prel = interfaceP[i]/meyerHardness;
            scalar km = 2*k[i]*nbrK[i]/(k[i] + nbrK[i]);
            scalar R = sqrt(sqr(R_[i]) + sqr(nbrR[i]));
            scalar Rf = max(0.5*(R_[i] + nbrR[i]), SMALL);
            scalar Rmult;

            if (Prel <= 0.0087)
            {
                Rmult = 333.3*Prel;
            }
            else
            {
                Rmult = 2.9 ;
            }

            scalar E = exp(5.738 - 0.528*log(3.937e7*Rf));

            if (Prel > 0.003)
            {
                hContact[i] = 0.4166*km*Prel*Rmult/R/E;
            }
            else if (Prel > 9e-6)
            {
                hContact[i] = 0.00125*km/R/E;
            }
            else
            {
                hContact[i] = 0.4166*km*sqrt(Prel)/R/E;
            }

        }
    }

    hGap_ = relax_*(hGas + hRad + hContact) + (1 - relax_)*hGap_;

    return (hGap_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::trisoGapFvPatchScalarField::
trisoGapFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    resistiveGapFvPatchScalarField(p, iF),
    gapWidth_(this->size(), 0),
    hGap_(this->size(), 0),
    R_(this->size(), 0),
    jd_(this->size(), -1),
    emissivityName_("emissivity"),
    gapGasModelName_("gapGas"),
    relax_(1),
    relaxGapWidth_(1),
    r1_(this->size(),0),
    r2_(this->size(),0)
{}


Foam::trisoGapFvPatchScalarField::
trisoGapFvPatchScalarField
(
    const trisoGapFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    resistiveGapFvPatchScalarField(ptf, p, iF, mapper),
    gapWidth_(ptf.gapWidth_),
    hGap_(ptf.hGap_),
    R_(ptf.R_),
    jd_(ptf.jd_),
    emissivityName_(ptf.emissivityName_),
    gapGasModelName_(ptf.gapGasModelName_),
    relax_(ptf.relax_),
    relaxGapWidth_(ptf.relaxGapWidth_)
{}


Foam::trisoGapFvPatchScalarField::
trisoGapFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    resistiveGapFvPatchScalarField(p, iF, dict, false),
    gapWidth_(this->size(), 0),
    hGap_(this->size(), 0),
    R_("roughness", dict, p.size()),
    jd_(this->size(), -1),
    emissivityName_(dict.lookupOrDefault<word>("emissivity", "emissivity")),
    gapGasModelName_(dict.lookupOrDefault<word>("gapGasModel", "gapGas")),
    relax_(dict.lookupOrDefault<scalar>("relax", 1.0)),
    relaxGapWidth_(dict.lookupOrDefault<scalar>("relaxGapWidth", 1.0))
{
    if( dict.found("alpha") )
    {
        alpha_gap_ = scalarField("alpha", dict, p.size());
    }

    if( dict.found("jumpDistance") )
    {
        jd_ = scalarField("jumpDistance", dict, p.size());
    }


}


Foam::trisoGapFvPatchScalarField::
trisoGapFvPatchScalarField
(
    const trisoGapFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    resistiveGapFvPatchScalarField(ptf, iF),
    gapWidth_(ptf.gapWidth_),
    hGap_(ptf.hGap_),
    R_(ptf.R_),
    jd_(ptf.jd_),
    emissivityName_(ptf.emissivityName_),
    gapGasModelName_(ptf.gapGasModelName_),
    relax_(ptf.relax_),
    relaxGapWidth_(ptf.relaxGapWidth_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::trisoGapFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    resistiveGapFvPatchScalarField::autoMap(m);
#ifdef OPENFOAMFOUNDATION    
    m(gapWidth_, gapWidth_);
    m(hGap_, hGap_);
    m(R_, R_);
    m(jd_, jd_);
#elif OPENFOAMESI
    gapWidth_.autoMap(m);
    hGap_.autoMap(m);
    R_.autoMap(m);
    jd_.autoMap(m);
#endif    
}


void Foam::trisoGapFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    resistiveGapFvPatchScalarField::rmap(ptf, addr);

    const trisoGapFvPatchScalarField& dmptf =
        refCast<const trisoGapFvPatchScalarField>(ptf);

    gapWidth_.rmap(dmptf.gapWidth_, addr);
    R_.rmap(dmptf.R_, addr);
    jd_.rmap(dmptf.jd_, addr);
}


void Foam::trisoGapFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const regionCoupledOFFBEATFvPatch& patch
        = refCast<const regionCoupledOFFBEATFvPatch>(regionCoupledPatch_);

    const regionCoupledOFFBEATFvPatch& nbrPatch
        = refCast<const regionCoupledOFFBEATFvPatch>(patch.nbrPatch());

    const trisoGapFvPatchScalarField& nbr = this->neighbour();

    if (regionCoupledPatch_.owner())
    {
        scalarField oldgapwidth = gapWidth_;
        // Update the gap width
        gapWidth_ =
        (
            0.5*
            (
                gapWidth()
              + patch.regionCoupledPatch().interpolate(nbr.gapWidth())
            )
        ) ;

        forAll(gapWidth_, faceI)
        {
            scalar r = 1;
            if(gapWidth_[faceI] > 5e-5)
            {
                r = 1;
            }
            else if (gapWidth_[faceI] > 0.5e-5)
            {
                r = -0.5e-5/4.5e-5 + relaxGapWidth_*5/4.5 + (1 - relaxGapWidth_)*gapWidth_[faceI]/(4.5e-5);
            }
            else
            {
                r = relaxGapWidth_;
            }

            gapWidth_[faceI] = r*gapWidth_[faceI] + (1-r)*oldgapwidth[faceI];

        }

        // Update the neighbour gap width
        nbr.gapWidth_ = nbrPatch.regionCoupledPatch().interpolate(gapWidth_);

        //- Update the interface heat transfer coefficient
        alpha() = 0.5*
        (
            hGap()
          + patch.regionCoupledPatch().interpolate(nbr.hGap())
        );
    }

    const fvPatchVectorField& totalDispPatch =
    (
        patch.boundaryMesh().mesh().foundObject<fvMesh>("referenceMesh")
    ) ?
    patch.lookupPatchField<volVectorField, vector>("DD")
    :
    patch.lookupPatchField<volVectorField, vector>("D");

    const fvPatchVectorField& totalDispNbrPatch =
    (
        patch.boundaryMesh().mesh().foundObject<fvMesh>("referenceMesh")
    ) ?
    nbrPatch.lookupPatchField<volVectorField, vector>("DD")
    :
    nbrPatch.lookupPatchField<volVectorField, vector>("D");

    vectorField nf = patch.Sf() / patch.magSf();

    vectorField Cf = patch.Cf()
                   + totalDispPatch;

    vectorField nbrCf = regionCoupledPatch_.regionCoupledPatch().interpolate
                      (
                             nbrPatch.Cf()
                           + totalDispNbrPatch
                      );

    resistiveGapFvPatchScalarField::updateCoeffs();
}


void Foam::trisoGapFvPatchScalarField::write(Ostream& os) const
{
    resistiveGapFvPatchScalarField::write(os);

#ifdef OPENFOAMFOUNDATION        
    writeEntryIfDifferent<word>(os, "emissivity", "emissivity", emissivityName_);
    writeEntryIfDifferent<word>(os, "gapGasModel", "gapGas", gapGasModelName_);
    writeEntry(os, "roughness", R_);
    writeEntry(os, "jumpDistance", jd_);
    writeEntryIfDifferent<scalar>(os, "relax", 1, relax_);
    writeEntryIfDifferent<scalar>(os, "relaxGapWidth", 1, relaxGapWidth_);
#elif OPENFOAMESI
    os.writeEntryIfDifferent<word>("emissivity", "emissivity", emissivityName_);
    os.writeEntryIfDifferent<word>("gapGasModel", "gapGas", gapGasModelName_);
    os.writeEntry("roughness", R_);
    os.writeEntry("jumpDistance", jd_);
    os.writeEntryIfDifferent<scalar>("relax", 1, relax_);
    os.writeEntryIfDifferent<scalar>("relaxGapWidth", 1, relaxGapWidth_);
#endif    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        trisoGapFvPatchScalarField
    );
};


// ************************************************************************* //
