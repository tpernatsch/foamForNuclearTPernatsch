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
#include "fuelRodGapFvPatchScalarField.H"
#include "regionCoupledOFFBEATFvPatch.H"
#include "gapGasModel.H"
#include "physicoChemicalConstants.H"
#include "gapContactFvPatchVectorField.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    const char* fuelRodGapFvPatchScalarField::group_ = 
        ("boundaryConditions::fuelRodGap");

    defineParameter(fuelRodGapFvPatchScalarField, F_hGas_, 
        "F_hGas", (dimless), 1.0);
    defineParameter(fuelRodGapFvPatchScalarField, F_hRad_, 
        "F_hRad", (dimless), 1.0);
    defineParameter(fuelRodGapFvPatchScalarField, F_hContact_, 
        "F_hContact", (dimless), 1.0);
        
    defineParameter(fuelRodGapFvPatchScalarField, delta_hGas_, 
        "delta_hGas", (dimless), 0.0);        
    defineParameter(fuelRodGapFvPatchScalarField, delta_hRad_, 
        "delta_hRad", (dimless), 0.0);        
    defineParameter(fuelRodGapFvPatchScalarField, delta_hContact_, 
        "delta_hContact", (dimless), 0.0);
}

// * * * * * * * * * * * * * * * * Private members  * * * * * * * * * * * * *//


Foam::tmp<Foam::scalarField> 
Foam::fuelRodGapFvPatchScalarField::gapWidth() const
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

    return max((nbrCf - Cf) & nf, 0.0);
}


Foam::tmp<Foam::scalarField> 
Foam::fuelRodGapFvPatchScalarField::hGap() const
{
    // Scaling factors
    const dimensionedScalar& F_hGas = this->F_hGas();
    const dimensionedScalar& F_hRad = this->F_hRad();
    const dimensionedScalar& F_hContact = this->F_hContact();
    const dimensionedScalar& delta_hGas = this->delta_hGas();
    const dimensionedScalar& delta_hRad = this->delta_hRad();
    const dimensionedScalar& delta_hContact = this->delta_hContact();

    const gapGasModel& gapGas = patch().boundaryMesh().mesh()
            .lookupObject<gapGasModel>(gapGasModelName_);
    
    const regionCoupledOFFBEATFvPatch& patch
        = refCast<const regionCoupledOFFBEATFvPatch>(regionCoupledPatch_);
    
    const regionCoupledOFFBEATFvPatch& nbrPatch
        = refCast<const regionCoupledOFFBEATFvPatch>(patch.nbrPatch());
    
    const fuelRodGapFvPatchScalarField& nbr = this->neighbour();
    
    const scalarField& pT = *this;
    scalarField nbrT = patch.regionCoupledPatch().interpolate(nbr);
    scalarField nbrR = patch.regionCoupledPatch().interpolate(nbr.R());
    
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
        {
            // Average interface temperature
            scalar T = 0.5*(pT[i] + nbrT[i]);            
            
            // scalar Tgas = gapGas.T();
            
            // Gap gas conductivity
            const scalar k = gapGas.kappa(T);
            
            // Interface pressure assumed zero (i.e. open gap)
            scalar pI = interfaceP[i]/1e4/9.8;

            // Gap width
            scalar d_gap = max(gapWidth_[i], 0);
            
            // roughness distance
            scalar d_eff = 0;
            
            d_eff = exp(-1.25e-3*pI)*(R_[i] + nbrR[i]);
            
            // Temperature jump distance
            // Info << "k: "<< k << ", T: "<<T << ", gapGas.p(): " << gapGas.p() << " ";
            scalar d_jump = 0.0137*k*sqrt(T) / (gapGas.p()*gapGas.a(T));

            scalar d_gap_total = max(d_gap + d_eff - 1.397e-6, 0.0 );
            
            hGas[i] = k / (d_gap_total + 1.8*d_jump);
            // hGas[i] = k / (d_gap + 1.8*d_jump + d_eff );
            hGas[i] = F_hGas.value()*max(hGas[i], 0.0) + delta_hGas.value();
            
        }
        // Radiative HTC
        {
            const scalar T1 = pT[i];
            const scalar T2 = nbrT[i];
            const scalar eps1 = max(eps[i], SMALL);
            const scalar eps2 = max(nbrEps[i], SMALL);
        
            hRad[i] = F_hRad.value()*(constant::physicoChemical::sigma.value()*
                      (T1 + T2)*(T1*T1 + T2*T2)/(1/eps1 + 1/eps2 - 1)) + delta_hRad.value();
        }
        // Contact HTC
        {
            const scalar meyerHardness = 680e6;
            
            scalar Prel = interfaceP[i]/meyerHardness;
            scalar km = 2*k[i]*nbrK[i]/(k[i] + nbrK[i]);
            scalar R = sqrt(sqr(R_[i]) + sqr(nbrR[i]));
            scalar Rf = 0.5*(R_[i] + nbrR[i]);
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

            hContact[i] = F_hContact.value()*hContact[i] + delta_hContact.value();
        }
    }
    
    hGap_ = relax_*(hGas + hRad + hContact) + (1 - relax_)*hGap_;

    return (hGap_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fuelRodGapFvPatchScalarField::
fuelRodGapFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    resistiveGapFvPatchScalarField(p, iF),
    gapWidth_(this->size(), 0),
    hGap_(this->size(), 0),
    R_(this->size(), 0),
    emissivityName_("emissivity"),
    gapGasModelName_("gapGas"),
    relax_(1)
{}


Foam::fuelRodGapFvPatchScalarField::
fuelRodGapFvPatchScalarField
(
    const fuelRodGapFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    resistiveGapFvPatchScalarField(ptf, p, iF, mapper),
    gapWidth_(mapper(ptf.gapWidth_)),
    hGap_(mapper(ptf.hGap_)),
    R_(mapper(ptf.R_)),
    emissivityName_(ptf.emissivityName_),
    gapGasModelName_(ptf.gapGasModelName_),
    relax_(ptf.relax_)
{}


Foam::fuelRodGapFvPatchScalarField::
fuelRodGapFvPatchScalarField
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
    emissivityName_(dict.lookupOrDefault<word>("emissivity", "emissivity")),
    gapGasModelName_(dict.lookupOrDefault<word>("gapGasModel", "gapGas")),
    relax_(dict.lookupOrDefault<scalar>("relax", 1.0))
{
    if( dict.found("alpha") )
    {
        alpha_gap_ = scalarField("alpha", dict, p.size());
    }
}


Foam::fuelRodGapFvPatchScalarField::
fuelRodGapFvPatchScalarField
(
    const fuelRodGapFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    resistiveGapFvPatchScalarField(ptf, iF),
    gapWidth_(ptf.gapWidth_),
    hGap_(ptf.hGap_),
    R_(ptf.R_),
    emissivityName_(ptf.emissivityName_),
    gapGasModelName_(ptf.gapGasModelName_),
    relax_(ptf.relax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fuelRodGapFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    resistiveGapFvPatchScalarField::autoMap(m);
#ifdef OPENFOAMFOUNDATION    
    m(gapWidth_, gapWidth_);
    m(hGap_, hGap_);
    m(R_, R_);
#elif OPENFOAMESI
    gapWidth_.autoMap(m);
    hGap_.autoMap(m);
    R_.autoMap(m);
#endif    
}


void Foam::fuelRodGapFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    resistiveGapFvPatchScalarField::rmap(ptf, addr);

    const fuelRodGapFvPatchScalarField& dmptf =
        refCast<const fuelRodGapFvPatchScalarField>(ptf);

    gapWidth_.rmap(dmptf.gapWidth_, addr);
    hGap_.rmap(dmptf.hGap_, addr);
    R_.rmap(dmptf.R_, addr);
}


void Foam::fuelRodGapFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const regionCoupledOFFBEATFvPatch& patch
        = refCast<const regionCoupledOFFBEATFvPatch>(regionCoupledPatch_);
    
    const regionCoupledOFFBEATFvPatch& nbrPatch
        = refCast<const regionCoupledOFFBEATFvPatch>(patch.nbrPatch());
    
    const fuelRodGapFvPatchScalarField& nbr = this->neighbour();
    
    if (regionCoupledPatch_.owner())
    {
        // Update the gap width
        gapWidth_ = 
        (
            0.5*
            (
                gapWidth()
              + patch.regionCoupledPatch().interpolate(nbr.gapWidth())
            )
        ) ;

        // Update the neighbour gap width
        nbr.gapWidth_ = nbrPatch.regionCoupledPatch().interpolate(gapWidth_);
        
        //- Update the interface heat transfer coefficient
        alpha() = 0.5*
        (
            hGap()
          + patch.regionCoupledPatch().interpolate(nbr.hGap())
        );
    }
    
    resistiveGapFvPatchScalarField::updateCoeffs();
}


void Foam::fuelRodGapFvPatchScalarField::write(Ostream& os) const
{
    resistiveGapFvPatchScalarField::write(os);
    
#ifdef OPENFOAMFOUNDATION        
    writeEntryIfDifferent<word>(os, "emissivity", "emissivity", emissivityName_);
    writeEntryIfDifferent<word>(os, "gapGasModel", "gapGas", gapGasModelName_);
    writeEntry(os, "roughness", R_);
    writeEntryIfDifferent<scalar>(os, "relax", 1, relax_);
#elif OPENFOAMESI
    os.writeEntryIfDifferent<word>("emissivity", "emissivity", emissivityName_);
    os.writeEntryIfDifferent<word>("gapGasModel", "gapGas", gapGasModelName_);
    R_.writeEntry("roughness", os);
    os.writeEntryIfDifferent<scalar>("relax", 1, relax_);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fuelRodGapFvPatchScalarField
    );
};


// ************************************************************************* //