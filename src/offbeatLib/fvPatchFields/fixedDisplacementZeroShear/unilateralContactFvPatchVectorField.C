/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "unilateralContactFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

unilateralContactFvPatchVectorField::
unilateralContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0),
    totalTraction_(p.size(), vector::zero),
    totalTraction0_(),
    relax_(1),
    // w_(p.size(), 1.0),
    currentTime_(db().time().value()),
    planeNormal_(vector::zero),
    planePoint_(vector::zero),
    contactRelRampWidth_(SMALL),
    contactRelOffset_(SMALL)
{}


unilateralContactFvPatchVectorField::
unilateralContactFvPatchVectorField
(
    const unilateralContactFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    solidDirectionMixedFvPatchVectorField(ptf, p, iF, mapper),
    traction_(ptf.traction_),
    pressure_(ptf.pressure_),
    totalTraction_(ptf.totalTraction_),
    totalTraction0_(ptf.totalTraction0_),
    relax_(ptf.relax_),
    // w_(ptf.w_),
    currentTime_(ptf.currentTime_),
    planeNormal_(ptf.planeNormal_),
    planePoint_(ptf.planePoint_),
    contactRelRampWidth_(ptf.contactRelRampWidth_),
    contactRelOffset_(ptf.contactRelOffset_),
    penaltyScaleFactor_(ptf.penaltyScaleFactor_)
{
    w_ = 0;

    if (mappingRequired)
    {
        if (mapper.hasUnmapped())
        {
            traction_ = vector::zero;
            totalTraction_ = vector::zero;
            pressure_ = 0.0;
        }

    #ifdef OPENFOAMFOUNDATION
        mapper(traction_, ptf.traction_);
        mapper(totalTraction_, ptf.totalTraction_);
        mapper(pressure_, ptf.pressure_);
    #elif OPENFOAMESI
        traction_.autoMap(mapper);
        totalTraction_.autoMap(mapper);
        pressure_.autoMap(mapper);
    #endif
    }
}


unilateralContactFvPatchVectorField::
unilateralContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0),
    totalTraction_(p.size(), vector::zero),
    totalTraction0_(),
    relax_(dict.lookupOrDefault<scalar>("relax", 1)),
    // w_(p.size(), 1.0),
    currentTime_(db().time().value()),
    planeNormal_(vector::zero),
    planePoint_(vector::zero),
    contactRelRampWidth_(dict.lookupOrDefault<scalar>(
        "contactRelRampWidth", 
        dict.lookupOrDefault<scalar>("relTolerance", SMALL)
    )),
    contactRelOffset_(dict.lookupOrDefault<scalar>("contactRelOffset", SMALL)),
    penaltyScaleFactor_(dict.lookupOrDefault<scalar>("penaltyScaleFactor", 1))
{
    //- Read refValue if present (otherwise equal)
    if (dict.found("refValue"))
    {
        refValue() =
        (
            vectorField("refValue", dict, p.size())
        );
    }
    else
    {
        refValue() =
        (
            vectorField("value", dict, p.size())
        );
    }

    //- Read gradient if present (otherwise zero gradient)
    if (dict.found("refGradient"))
    {
        refGrad() =
        (
            vectorField("refGradient", dict, p.size())
        );
    }
    else
    {
        refGrad() = vector::zero;
    }

    if (dict.found("traction"))
    {
        traction_ =         
        (
            vectorField("traction", dict, p.size())
        );
    }

    if (dict.found("pressure"))
    {
        pressure_ =         
        (
            scalarField("pressure", dict, p.size())
        );
    }

    if(patch().boundaryMesh().mesh().foundObject<gapGasModel>("gapGas"))
    {
        const gapGasModel& gapGas
        = patch().boundaryMesh().mesh().lookupObject<gapGasModel>("gapGas");

        pressure_ = gapGas.p();
    }
    
    // Read totalTraction field if present in dictionary
    if (dict.found("totalTraction"))
    {
        totalTraction_ = vectorField("totalTraction" , dict, p.size());
    }

    if (dict.found("planeNormal"))
    {
#ifdef OPENFOAMFOUNDATION        
        planeNormal_ = dict.lookup<vector>("planeNormal");
#elif OPENFOAMESI
        planeNormal_ = dict.get<vector>("planeNormal");
#endif    
    }
    else
    {
        planeNormal_ = gSum(patch().nf()*patch().magSf())/gSum(patch().magSf());
    }

    if (dict.found("planePoint"))
    {
        planePoint_ = dict.lookupOrDefault<vector>("planePoint", vector::zero);
    }
    else
    {
        if(p.size())
        {     
            vectorField Cf = patch().Cf();
            
            if (db().foundObject<fvMesh>("referenceMesh"))
            {
                const fvMesh& refMesh =  db().lookupObject<fvMesh>("referenceMesh");
                Cf = refMesh.C().boundaryField()[patch().index()];
            }
            
            scalar minDist = (Cf[0] & (-planeNormal_));
            int index=0;
            
            forAll(Cf, faceI)
            {
                if( (Cf[faceI] & (-planeNormal_)) < minDist)
                {
                    minDist = (Cf[faceI] & (-planeNormal_));
                    index = faceI;
                }
            }
            
            planePoint_ = Cf[index];
        }
    }
    
    // Crate old totalTraction vector for DD solvers
    if(internalField().name() == "DD")
    {
#ifdef OPENFOAMFOUNDATION         
        totalTraction0_.set(new vectorField(p.size(), vector::zero));
#elif OPENFOAMESI
        totalTraction0_.reset(new vectorField("totalTraction", dict, p.size()));
#endif      
    }

    this->valueFraction() = sqr(planeNormal_);
    
    Field<vector> normalValue = transform(valueFraction(), refValue());

    Field<vector> gradValue =
        this->patchInternalField() + refGrad()/this->patch().deltaCoeffs();

    const symmTensorField valueFractionCopy = valueFraction();
    Field<vector> transformGradValue = 
        transform(I - valueFractionCopy, gradValue);

    Field<vector>::operator=(normalValue + transformGradValue);

    if(dict.found("weights"))
    {
        w_ = scalarField("weights" , dict, p.size());
    }
    // else
    // {
    //     // Penetration with sign (positive if penetrating)
    //     scalarField pen = (Df + Cf - planePoint_) & planeNormal_;   

    //     // Cap penetration pen >= 0
    //     pen = max(pen, scalar(0.0));

    //     // Smoothing distance
    //     scalarField deltaPen = contactRelRampWidth_/patch().deltaCoeffs();
        
    //     // Smooth coefficient, goes from 0 (just touching) to 1 (penetration >= deltaPen)
    //     scalarField s = pen / deltaPen;
        
    //     // clamp: s ∈ [0,1]
    //     s = min(s, scalar(1.0));  

        
    //     // s = 0          -> wNew = 1.0 (no contact)
    //     // s = 1          -> wNew = contactWeight (full contact)
    //     // 0 < s < 1      -> smooth transition inside contact
    //     scalarField contactWeight = 1.0 - Rp/(Rp + Rs);
    //     scalarField wNew = (1.0 - s)*1.0 + s*contactWeight;
        
    //     // Relaxation
    //     w_ = relax_*wNew + (1.0 - relax_)*w_;
    // }
}


unilateralContactFvPatchVectorField::
unilateralContactFvPatchVectorField
(
    const unilateralContactFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(ptf, iF),
    traction_(ptf.traction_),
    pressure_(ptf.pressure_),
    totalTraction_(ptf.totalTraction_),
    totalTraction0_(),
    relax_(ptf.relax_),
    // w_(ptf.w_),
    currentTime_(ptf.currentTime_),
    planeNormal_(ptf.planeNormal_),
    planePoint_(ptf.planePoint_),
    contactRelRampWidth_(ptf.contactRelRampWidth_),
    contactRelOffset_(ptf.contactRelOffset_),
    penaltyScaleFactor_(ptf.penaltyScaleFactor_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void unilateralContactFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidDirectionMixedFvPatchVectorField::autoMap(m);
#ifdef OPENFOAMFOUNDATION
    m(traction_, traction_);
    m(pressure_, pressure_);
    m(totalTraction_, totalTraction_);
#elif OPENFOAMESI    
    traction_.autoMap(m);
    totalTraction_.autoMap(m);
    pressure_.autoMap(m);
#endif    
}


// Reverse-map the given fvPatchField onto this fvPatchField
void unilateralContactFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    solidDirectionMixedFvPatchVectorField::rmap(ptf, addr);

    const unilateralContactFvPatchVectorField& dmptf =
        refCast<const unilateralContactFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
    totalTraction_.rmap(dmptf.totalTraction_, addr);
}


void unilateralContactFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    
    // Find whether or not the simulation usese a small-strain incremental
    // solver with mesh update, to apply a special treatement
    bool smallStrainIncUpdated = false;
    if( (internalField().name() == "DD") and 
        not(patch().boundaryMesh().mesh().foundObject<volTensorField>("Finv")))
    {
        smallStrainIncUpdated = true;
    }
    
    // Update old totalTraction vector
    if(db().time().value() > currentTime_)
    {
        if(smallStrainIncUpdated)
        {
            totalTraction0_() = totalTraction_;
        }
        currentTime_ = db().time().value();
    }

    // Patch normal 
    vectorField n(patch().nf());
    vectorField nCurrent(n);
    
    if(patch().boundaryMesh().mesh().foundObject<volTensorField>("Finv"))
    {
        const tensorField& Finv =
        patch().boundaryMesh().mesh().foundObject<fvMesh>("referenceMesh") ?
        patch().lookupPatchField<volTensorField, tensor>("relFinv") :
        patch().lookupPatchField<volTensorField, tensor>("Finv");
        
        // Patch unit normals (deformed configuration)
        nCurrent = (Finv.T() & n);
        nCurrent /= mag(nCurrent);
    };


    // // Include the gap gas pressure
    // if(patch().boundaryMesh().mesh().foundObject<gapGasModel>("gapGas"))
    // {
    //     const gapGasModel& gapGas
    //     = patch().boundaryMesh().mesh().lookupObject<gapGasModel>("gapGas");

    //     pressure_ = gapGas.p();
    // }

    totalTraction_ = traction_ - pressure_*nCurrent;
    vectorField deltaTraction(patch().size(), vector::zero);
    if(smallStrainIncUpdated)
    {        
        deltaTraction = totalTraction_ - totalTraction0_();   

        // Absolute threshold in Pascals (adjust based on your case)
        scalar absTol = 1e-6;  
        scalar relTol = 10 * SMALL;  

        forAll(deltaTraction, i)
        {
            scalar refValue = mag(totalTraction_[i]);
            scalar tol = max(absTol, relTol * refValue); 

            if (mag(deltaTraction[i]) < tol)
            {
                deltaTraction[i] = vector::zero;
            }
        }
    }

    // Set gradient to force zero shear traction
    const fvPatchField<scalar>& mu =
        patch().lookupPatchField<volScalarField, scalar>("mu");

    const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>("lambda");
        
    scalarField K =
            patch().lookupPatchField<volScalarField, scalar>("threeK")/3;
        
    const fvPatchField<tensor>& sigmaExp = 
        smallStrainIncUpdated ?
        patch().lookupPatchField<volTensorField, tensor>("DSigmaExp") :
        patch().lookupPatchField<volTensorField, tensor>("sigmaExp");
    
    // Set 2mu+lambda field
    scalarField twoMuLambda(2*mu + lambda);
    
    // Patch resistance
    scalarField Rp = 1/patch().deltaCoeffs()/max(twoMuLambda, SMALL);

    // Spring-contact resistance (penalty-like)
    scalarField Rs = 1/patch().deltaCoeffs()/max(penaltyScaleFactor_*K, SMALL);

    // Calculate normal component of displacement
    vectorField Df = smallStrainIncUpdated ?
    (patchInternalField() - Rp*(n & sigmaExp) + Rp*deltaTraction) :
    (patchInternalField() - Rp*(n & sigmaExp) + Rp*totalTraction_);
    
    if (internalField().name() == "DD")
    {
        const fvPatchVectorField& Doldf =
        db().lookupObject<volVectorField>("D").oldTime().boundaryField()[patch().index()];

        Df += Doldf;
    }

    vectorField Cf =
    db().foundObject<fvMesh>("referenceMesh") ?
    db().lookupObject<fvMesh>("referenceMesh").boundaryMesh()[patch().index()].faceCentres() :
    patch().Cf();

    // Physical signed gap
    scalarField gap = ((Df + Cf - planePoint_) & planeNormal_);

    // activation offset (dimensionless -> length via deltaCoeffs)
    scalarField offset = contactRelOffset_/patch().deltaCoeffs();
    
    // effective penetration-like measure (with positive/negative offset)
    scalarField pen = max(gap + offset, scalar(0.0));

    // Smoothing distance
    scalarField deltaPen(contactRelRampWidth_/(patch().deltaCoeffs()));
    
    // Smooth coefficient, goes from 0 (just touching) to 1 (penetration >= deltaPen)
    scalarField s = pen / deltaPen;
    
    // clamp: s ∈ [0,1]
    s = min(s, scalar(1.0));  
    
    // s = 0          -> wNew = 1.0 (no contact)
    // s = 1          -> wNew = contactWeight (full contact)
    // 0 < s < 1      -> smooth transition inside contact
    scalarField contactWeight = 1.0 - Rp/(Rp + Rs);
    scalarField wNew = (1.0 - s)*1.0 + s*contactWeight;
    
    // Relaxation
    w_ = relax_*wNew + (1.0 - relax_)*w_;

    // Set value for normal component
    // NOTE: For w_>0 the patch must be equivalent to fixed gradient
    // thus the internalField is added at the evaluate and snGrad stage
    if(internalField().name() == "DD")
    {
        const fvPatchVectorField& Doldf =
        db().lookupObject<volVectorField>("D").oldTime().boundaryField()[patch().index()];
        // DD = incremental disp. to plane if w_ = 0
        if (smallStrainIncUpdated)
        {
           refValue() = (1 - w_)*(planePoint_ - Cf - Doldf) + w_*(-Rp*(n & sigmaExp) + Rp*deltaTraction);
        }
        else
        {
            refValue() = (1 - w_)*(planePoint_ - Cf - Doldf) + w_*(-Rp*(n & sigmaExp) + Rp*totalTraction_);
        }
    }
    else
    {
        // DD = incremental disp to plane so D = planePoint_ - Cf if w_ = 0
        refValue() = (1 - w_)*(planePoint_ - Cf) + (w_)*(-Rp*(n & sigmaExp) + Rp*totalTraction_);
    }

    // Set gradient for tangent component
    refGrad() = smallStrainIncUpdated ? 
    (w_*deltaTraction - (n & sigmaExp)) / twoMuLambda :
    (w_*totalTraction_ - (n & sigmaExp)) / twoMuLambda;

    solidDirectionMixedFvPatchVectorField::updateCoeffs();
}


// Write
void unilateralContactFvPatchVectorField::write(Ostream& os) const
{
    solidDirectionMixedFvPatchVectorField::write(os);
    
#ifdef OPENFOAMFOUNDATION  
    writeEntry<scalarField>(os, "pressure", pressure_);
#elif OPENFOAMESI
    pressure_.writeEntry( "pressure", os);
#endif

#ifdef OPENFOAMFOUNDATION  
    writeEntry<vectorField>(os, "traction", traction_);
#elif OPENFOAMESI
    traction_.writeEntry( "traction", os);
#endif 


#ifdef OPENFOAMFOUNDATION  
    writeEntry<vectorField>(os, "totalTraction", totalTraction_);
#elif OPENFOAMESI
    totalTraction_.writeEntry( "totalTraction", os);
#endif 

#ifdef OPENFOAMFOUNDATION 

    writeEntryIfDifferent<scalar>(os, "relax", 1, relax_);
    writeEntryIfDifferent<scalar>(os, "contactRelRampWidth", SMALL, contactRelRampWidth_);
    writeEntryIfDifferent<scalar>(os, "contactRelOffset", SMALL, contactRelOffset_);
    writeEntry<vector>(os, "planePoint", planePoint_);
    writeEntry<vector>(os, "planeNormal", planeNormal_);
    writeEntry(os, "weights", w_);

#elif OPENFOAMESI

    os.writeEntryIfDifferent<scalar>("relax", 1, relax_);
    os.writeEntryIfDifferent<scalar>("contactRelRampWidth", SMALL, contactRelRampWidth_);
    os.writeEntryIfDifferent<scalar>("contactRelOffset", SMALL, contactRelOffset_);
    os.writeEntry("planePoint", planePoint_);
    os.writeEntry("planeNormal", planeNormal_);
    w_.writeEntry("weights", os);

#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    unilateralContactFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
