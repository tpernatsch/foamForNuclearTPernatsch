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

#include "tractionDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    totalTraction_(p.size(), vector::zero),
    totalTraction0_(),
    relax_(1),
    stressName_("sigma"),
    planeStrain_(false),
    flatSurface_(false),
    pressureList_(),
    tractionList_(),
    fixedSpring_(false),
    fixedSpringModulus_(3.5e3),
    dashpotModulus_(3.5e3),
    // w_(p.size(), 1.0),
    currentTime_(db().time().value())
{}


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    solidDirectionMixedFvPatchVectorField(tdpvf, p, iF, mapper),
    traction_(tdpvf.traction_, mapper),
    pressure_(tdpvf.pressure_, mapper),
    totalTraction_(tdpvf.totalTraction_, mapper),
    totalTraction0_(),
    relax_(tdpvf.relax_),
    stressName_(tdpvf.stressName_),
    planeStrain_(tdpvf.planeStrain_),
    flatSurface_(tdpvf.flatSurface_),
    pressureList_(tdpvf.pressureList_),
    tractionList_(tdpvf.tractionList_),
    fixedSpring_(tdpvf.fixedSpring_),
    fixedSpringModulus_(tdpvf.fixedSpringModulus_),
    dashpotModulus_(tdpvf.dashpotModulus_),
    // w_(tdpvf.w_),
    currentTime_(tdpvf.currentTime_)
{
    if (mappingRequired)
    {
        if (mapper.hasUnmapped())
        {
            traction_ = vector::zero;
            pressure_ = 0.0;
            totalTraction_ = vector::zero;
        }

    #ifdef OPENFOAMFOUNDATION
        mapper(traction_, tdpvf.traction_);
        mapper(pressure_, tdpvf.pressure_);
        mapper(totalTraction_, tdpvf.totalTraction_);
    #endif
    }
 }


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
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
    stressName_(dict.lookupOrDefault<word>("stressName", "sigma")),
    planeStrain_(dict.lookupOrDefault<bool>("planeStrain", false)),
    flatSurface_(dict.lookupOrDefault<bool>("flatSurface", false)),
    pressureList_(),
    tractionList_(),
    fixedSpring_(dict.lookupOrDefault<bool>("fixedSpring", false)),
    fixedSpringModulus_(fixedSpring_? readScalar(
        dict.lookup("fixedSpringModulus")) : 0.0),
    dashpotModulus_(fixedSpring_ ? readScalar(
        dict.lookup("dashpotModulus")) : 0.0),
    // w_(p.size(), 1.0),
    currentTime_(db().time().value())
{
    if(valueRequired)
    {
        // Read traction (either list or fixed traction)
        if(dict.found("tractionList"))
        {
#ifdef OPENFOAMFOUNDATION
            tractionList_.set
            (
                new Function1s::Table<vector>
                (
                    "tractionList", dict.subDict("tractionList")
                )
            );
#elif OPENFOAMESI
            word listCoeffsName = "tractionList";
            if(dict.found("tractionListCoeffs"))
            {
                listCoeffsName = "tractionListCoeffs";
            }
            tractionList_.reset
            (
                new vectorTable
                (
                    "tractionList", dict.subDict(listCoeffsName)
                )
            );
#endif

            traction_ =
            tractionList_->value
            (
                patch().boundaryMesh().mesh().time().timeToUserTime
                (
                    this->db().time().value()
                )
            );
        }
        else if (dict.found("traction"))
        {
            traction_ = vectorField("traction", dict, patch().size());
        }
        else
        {
            FatalErrorInFunction()
            << "tractionDisplacement BC for " << patch().name() << " requires:" << nl
            << "- \"tractionList\" or" << nl
            << "- \"traction\"" << abort(FatalError) << endl;
        }

        // Read pressure (either list or fixed pressure)
        if(dict.found("pressureList"))
        {
#ifdef OPENFOAMFOUNDATION
            pressureList_.set
            (
                new Function1s::Table<scalar>
                (
                    "pressureList", dict.subDict("pressureList")
                )
            );
#elif OPENFOAMESI
            word listCoeffsName = "pressureList";
            if(dict.found("pressureListCoeffs"))
            {
                listCoeffsName = "pressureListCoeffs";
            }
            pressureList_.reset
            (
                new scalarTable
                (
                    "pressureList", dict.subDict("pressureList")
                )
            );
#endif

            pressure_ =
            pressureList_->value
            (
                patch().boundaryMesh().mesh().time().timeToUserTime
                (
                    this->db().time().value()
                )
            );
        }
        else if (dict.found("pressure"))
        {
            pressure_ = scalarField("pressure", dict, patch().size());
        }
        else
        {
            FatalErrorInFunction()
            << "tractionDisplacement BC for " << patch().name() << " requires:" << nl
            << "- \"pressureList\" or" << nl
            << "- \"pressure\"" << abort(FatalError) << endl;
        }
    }

    // Read initial value
    refValue() = vectorField("value", dict, p.size());

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

    if (planeStrain_ and flatSurface_)
    {
        FatalErrorInFunction() 
        << "You cannot activate both \"planeStrain\" and \"flatSurface\" option for tractionDisplacement BC for " << patch().name() << nl
        << "Select only one of the two options." << nl
        << abort(FatalError) << endl;
    }
    else if (planeStrain_)
    {
        // Act as normal fixedGradient BC
        this->valueFraction() = symmTensor::zero;
        this->refValue() = vector::zero;
    }
    else if (flatSurface_ or fixedSpring_)
    {
        // Non-plane strain BC, with flat surface and/or fixed spring act as 
        // mixedDirection BC
        this->valueFraction() = sqr(patch().nf());
        if (flatSurface_)
        {
            w_ = 0.0;
        }
    }
    else
    {
        // Act as standard fixed gradient BC
        this->valueFraction() = symmTensor::zero;
        this->refValue() = vector::zero;
    }

    Field<vector> normalValue = transform(valueFraction(), refValue());

    Field<vector> gradValue =
        this->patchInternalField() + refGrad()/this->patch().deltaCoeffs();

    const symmTensorField valueFractionCopy = valueFraction();
    Field<vector> transformGradValue = 
        transform(I - valueFractionCopy, gradValue);

    Field<vector>::operator=(normalValue + transformGradValue);

    // Read spring-dashpot weights if applicable
    if(fixedSpring_ & dict.found("weights"))
    {
        w_ =
        (
            scalarField("weights", dict, p.size())
        );
    }    
    
    // Read totalTraction field if present in dictionary
    if (dict.found("totalTraction"))
    {
        totalTraction_ = vectorField("totalTraction" , dict, p.size());
    }
}


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:  
    solidDirectionMixedFvPatchVectorField(tdpvf, iF),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_),
    totalTraction_(tdpvf.totalTraction_),
    totalTraction0_(),
    relax_(tdpvf.relax_),
    stressName_(tdpvf.stressName_),
    planeStrain_(tdpvf.planeStrain_),
    flatSurface_(tdpvf.flatSurface_),
    pressureList_(),
    tractionList_(),
    fixedSpring_(false),
    fixedSpringModulus_(tdpvf.fixedSpringModulus_),
    dashpotModulus_(tdpvf.dashpotModulus_),
    // w_(tdpvf.w_),
    currentTime_(tdpvf.currentTime_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tractionDisplacementFvPatchVectorField::autoMap
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
    pressure_.autoMap(m);
    totalTraction_.autoMap(m);
#endif
}


void tractionDisplacementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    solidDirectionMixedFvPatchVectorField::rmap(ptf, addr);

    const tractionDisplacementFvPatchVectorField& dmptf =
        refCast<const tractionDisplacementFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
    totalTraction_.rmap(dmptf.totalTraction_, addr);
}


// Foam::tmp<Foam::Field<vector>> 
// Foam::tractionDisplacementFvPatchVectorField::snGrad() const
// {
//     // // Lookup the gradient field (D or DD in dependance of the solver used)
//     // const fvPatchField<tensor>& gradField =
//     //     patch().lookupPatchField<volTensorField, tensor>
//     //     (
//     //         "grad" + internalField().name()
//     //     );  

//     // // Face unit normal
//     // const vectorField n = patch().nf();

//     // // Cell centers and delta vectors
//     // vectorField ownCp(patch().boundaryMesh().mesh().C(), patch().faceCells());
//     // const vectorField delta = patch().Cf() - ownCp;

//     // // Non-orthogonal correction vectors
//     // vectorField k = delta - n*(n&delta);

//     // return gradient() + patch().deltaCoeffs()*(w_ - 1)*(
//     //     patchInternalField() + (k & gradField));

//     return gradient();
// }


void tractionDisplacementFvPatchVectorField::updateTraction()
{
    if(pressureList_.valid())
    {
        pressure_ =
        pressureList_->value
        (
            patch().boundaryMesh().mesh().time().timeToUserTime
            (
                this->db().time().value()
            )
        );
    }

    if(tractionList_.valid())
    {
        traction_ =
        tractionList_->value
        (
            patch().boundaryMesh().mesh().time().timeToUserTime
            (
                this->db().time().value()
            )
        );
    }
}


void tractionDisplacementFvPatchVectorField::updateCoeffs()
{
    if (updated())
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

    // Update traction
    updateTraction();

    // Update old totalTraction vector
    if(db().time().value() > currentTime_)
    {
        if(smallStrainIncUpdated)
        {
            // Crate old totalTraction vector for DD solvers
            if(!totalTraction0_.valid())
            {
                totalTraction0_.set(new vectorField(patch().size(), vector::zero));
            }
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

    // Face area and total area
    scalarField magSf(patch().magSf());
    scalar totalArea(gSum(magSf));

    // Patch number of faces
    scalar nFaces(patch().size());
    scalar nFacesTotal(nFaces);
    reduce(nFacesTotal, sumOp<scalar>());

    // Patch fields
    const fvPatchField<scalar>& mu =
        patch().lookupPatchField<volScalarField, scalar>("mu");

    const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>("lambda");

    scalarField twoMuLambda(2*mu + lambda);

    // Lookup the explicit component of the stress
    const fvPatchField<tensor>& sigmaExp = smallStrainIncUpdated ?
    patch().lookupPatchField<volTensorField, tensor>("DSigmaExp") :
    patch().lookupPatchField<volTensorField, tensor>("sigmaExp");

    // Local resistance (needed to simplify notation in planeStrain and
    // fixedSpring cases, but kept in all cases for consistency)
    scalarField Rp = 1/patch().deltaCoeffs()/(max(2*mu + lambda, VSMALL));

    // Boundary-specified force per unit area
    vectorField Fb = traction_ - pressure_*nCurrent;

    // Update total force on patch
    Ftot_ = gSum(Fb*magSf);

    // Calculate the new gradient R_P*deltaCoeff is equivalent to 1/twoMuLambda
    // for cases with no planeStrain or cases where RP is uniform
    vectorField newGradient(smallStrainIncUpdated ?
        Rp*patch().deltaCoeffs()*(deltaTraction - (n & sigmaExp)) :
        Rp*patch().deltaCoeffs()*(totalTraction_ - (n & sigmaExp)));

    // Force the normal component of the boundary displacement 
    // or the one of the normal boundary gradient to be uniform
    // TODO: write derivation in manual
    if(flatSurface_ or planeStrain_)
    {
        // Average 1/Rp
        scalar RpInvAvg = 1/totalArea*gSum(1/Rp * magSf);

        // Average field Dp/Rp
        vector DpOverRpAvg = 1/totalArea*gSum(patchInternalField()/Rp * magSf);
        
        // Average stress explicit component
        vector sigmaExpAvg(1/totalArea*gSum(n & sigmaExp * magSf));

        // Average external traction
        vector tractionAvg = smallStrainIncUpdated ?
            1/totalArea*gSum(deltaTraction * magSf) :
            1/totalArea*gSum(totalTraction_ * magSf);

        if (fixedSpring_)
        {
            // Spring & dashpot resistances
            scalarField Rs = magSf/max(fixedSpringModulus_/nFacesTotal, SMALL);
            scalarField Rd = magSf/max(dashpotModulus_/nFacesTotal, SMALL); 
            
            // Adjust average 1/Rp (i.e. 1/Rp + 1/Rs + 1/Rd)
            RpInvAvg += 1/totalArea*gSum((1/Rs + 1/Rd) * magSf);
            
            // Add dashpot traction to average traction
            tractionAvg += 1/totalArea*gSum((*this/Rd) * magSf);            
        }      
        
        newGradient = (newGradient & (I - sqr(n)));

        if(flatSurface_)
        {    
            // Uniform boundary normal displacement, mixedDirection BC
            refValue() = (1 - relax_)*refValue() + relax_*(1/RpInvAvg*(sqr(n) & (tractionAvg - sigmaExpAvg + DpOverRpAvg)));
        }
        else // plane strain
        {
            // Uniform boundary normal gradient (normal strain component), fixedGradient BC
            newGradient += patch().deltaCoeffs()*(1/RpInvAvg*(sqr(n) & (tractionAvg - sigmaExpAvg)));
        }
        
    }
    // No plane strain, nor flat surface. Only resistance network for spring-dashpot system
    else if(fixedSpring_)
    {
        // Spring resistance
        scalarField Rs = magSf/max(fixedSpringModulus_/nFacesTotal, SMALL);

        // Dashpot resistance
        scalarField Rd = magSf/max(dashpotModulus_/nFacesTotal, SMALL); 

        // Spring-dashpot resistance
        scalarField Rsd = Rs*Rd/(Rs + Rd);

        // Total resistance
        scalarField Rtot = Rp + Rsd;
        
        // Update weight factor
        w_ = relax_*(Rsd/Rtot) + (1 - relax_)*w_;

        // Fix normal value
        // NOTE: For w_>0 the patch must be equivalent to fixed gradient
        // thus the internalField is added at the evaluate and snGrad stage
        refValue() = smallStrainIncUpdated ?
            w_*(Rp*(deltaTraction - (n & sigmaExp) + (*this/Rd))) :
            w_*(Rp*(totalTraction_ - (n & sigmaExp) + (*this/Rd)));
    }
    else
    {
        refValue() = vector::zero;
        valueFraction() = symmTensor::zero;
        w_ = 1.0;
    }


    // Relax gradient
    refGrad() = (1-relax_)*refGrad() + relax_*newGradient;

    solidDirectionMixedFvPatchVectorField::updateCoeffs();
}

void tractionDisplacementFvPatchVectorField::write(Ostream& os) const
{
    solidDirectionMixedFvPatchVectorField::write(os);

    if (tractionList_.valid())
    {
#ifdef OPENFOAMFOUNDATION
        tractionList_->write(os);
#elif OPENFOAMESI
        tractionList_->writeData(os);
#endif
    }
    else
    {
#ifdef OPENFOAMFOUNDATION
        writeEntry<vectorField>(os, "traction", traction_);
#elif OPENFOAMESI
        traction_.writeEntry( "traction", os);
#endif
    }

    if (pressureList_.valid())
    {
#ifdef OPENFOAMFOUNDATION
        pressureList_->write(os);
#elif OPENFOAMESI
        pressureList_->writeData(os);
#endif
    }
    else
    {
#ifdef OPENFOAMFOUNDATION
        writeEntry<scalarField>(os, "pressure", pressure_);
#elif OPENFOAMESI
        pressure_.writeEntry( "pressure", os);
#endif
    }

#ifdef OPENFOAMFOUNDATION  
    writeEntry<vectorField>(os, "totalTraction", totalTraction_);
#elif OPENFOAMESI
    totalTraction_.writeEntry( "totalTraction", os);
#endif  

#ifdef OPENFOAMFOUNDATION    
    writeEntryIfDifferent<scalar>(os, "relax", 1, relax_);
    writeEntryIfDifferent<word>(os, "stressName", "sigma", stressName_);
    writeEntryIfDifferent<bool>(os, "planeStrain", false, planeStrain_);
    writeEntryIfDifferent<bool>(os, "flatSurface", false, flatSurface_);

    writeEntryIfDifferent<bool>(os, "fixedSpring", false, fixedSpring_);
    if(fixedSpring_)
    {
        writeEntry<scalar>(os, "fixedSpringModulus", fixedSpringModulus_);
        writeEntry<scalar>(os, "dashpotModulus", dashpotModulus_);
        writeEntry(os, "weights", w_);
    }

    writeEntry(os, "value", *this);
#elif OPENFOAMESI
    os.writeEntryIfDifferent<scalar>("relax", 1, relax_);
    os.writeEntryIfDifferent<word>("stressName", "sigma", stressName_);
    os.writeEntryIfDifferent<bool>("planeStrain", false, planeStrain_);
    os.writeEntryIfDifferent<bool>("flatSurface", false, flatSurface_);

    os.writeEntryIfDifferent<bool>("fixedSpring", false, fixedSpring_);
    if(fixedSpring_)
    {
        os.writeEntry<scalar>("fixedSpringModulus", fixedSpringModulus_);
        os.writeEntry<scalar>("dashpotModulus", dashpotModulus_);
        w_.writeEntry("weights", os);
    }


    this->writeEntry("value", os);
#endif

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    tractionDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
