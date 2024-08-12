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
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    relax_(1),
    stressName_("sigma"),
    planeStrain_(false),
    pressureList_(),
    tractionList_(),
    fixedSpring_(false),
    fixedSpringModulus_(3.5e3),
    dashpotModulus_(3.5e3),
    w_(p.size(), 1.0)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_),
    relax_(tdpvf.relax_),
    stressName_(tdpvf.stressName_),
    planeStrain_(tdpvf.planeStrain_),
    pressureList_(),
    tractionList_(),
    fixedSpring_(false),
    fixedSpringModulus_(tdpvf.fixedSpringModulus_),
    dashpotModulus_(tdpvf.dashpotModulus_),
    w_(tdpvf.w_)
{
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
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0),
    relax_(dict.lookupOrDefault<scalar>("relax", 1)),
    stressName_(dict.lookupOrDefault<word>("stressName", "sigma")),
    planeStrain_(dict.lookupOrDefault<bool>("planeStrain", false)),
    pressureList_(),
    tractionList_(),
    fixedSpring_(dict.lookupOrDefault<bool>("fixedSpring", false)),
    fixedSpringModulus_(fixedSpring_? readScalar(
        dict.lookup("fixedSpringModulus")) : 0.0),
    dashpotModulus_(fixedSpring_ ? readScalar(
        dict.lookup("dashpotModulus")) : 0.0),
    w_(p.size(), 1.0)
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
            tractionList_.reset
            ( 
                new vectorTable
                (
                    "tractionList", dict.subDict("tractionList")
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

    // Read gradient if present (otherwise zero gradient)
    if (dict.found("gradient"))
    {
        gradient() =
        (
            vectorField("gradient", dict, p.size())
        );
    }
    else
    {
        gradient() = vector::zero;
    }

    // Read displacement value if present (otherwise zero gradient)
    if (dict.found("value"))
    {
        fvPatchVectorField::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }

    // Read spring-dashpot weights if applicable
    if(fixedSpring_ & dict.found("weights"))
    {
        w_ =
        (
            scalarField("weights", dict, p.size())
        );
    }
}    


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:  
    fixedGradientFvPatchVectorField(tdpvf, iF),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_),
    relax_(tdpvf.relax_),
    stressName_(tdpvf.stressName_),
    planeStrain_(tdpvf.planeStrain_),
    pressureList_(),
    tractionList_(),
    fixedSpring_(false),
    fixedSpringModulus_(tdpvf.fixedSpringModulus_),
    dashpotModulus_(tdpvf.dashpotModulus_),
    w_(tdpvf.w_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tractionDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
#ifdef OPENFOAMFOUNDATION
    m(traction_, traction_);
    m(pressure_, pressure_);
#elif OPENFOAMESI    
    traction_.autoMap(m);
    pressure_.autoMap(m);
#endif    
}


void tractionDisplacementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const tractionDisplacementFvPatchVectorField& dmptf =
        refCast<const tractionDisplacementFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


Foam::tmp<Foam::Field<vector>> 
Foam::tractionDisplacementFvPatchVectorField::snGrad() const
{
    // Lookup the gradient field (D or DD in dependance of the solver used)
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad" + internalField().name()
        );  

    // Face unit normal
    const vectorField n = patch().nf();

    // Cell centers and delta vectors
    vectorField ownCp(patch().boundaryMesh().mesh().C(), patch().faceCells());
    const vectorField delta = patch().Cf() - ownCp;

    // Non-orthogonal correction vectors
    vectorField k = delta - n*(n&delta);

    return gradient() + patch().deltaCoeffs()*(w_ - 1)*(
        patchInternalField() + (k & gradField));
}


void tractionDisplacementFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Patch normal 
    vectorField n(patch().nf());
    vectorField nCurrent(n);

    // Face area and total area
    scalarField magSf(patch().magSf());
    scalar totalArea(gSum(magSf));

    // Patch number of faces
    scalar nFaces(patch().size());
    scalar nFacesTotal(nFaces);
    reduce(nFacesTotal, sumOp<scalar>());

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


    // Patch fields
    const fvPatchField<scalar>& mu =
        patch().lookupPatchField<volScalarField, scalar>("mu");

    const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>("lambda");
    
    scalarField twoMuLambda(2*mu + lambda);

    // Lookup the explicit component of the stress
    const fvPatchField<tensor>& sigmaExp =
        patch().lookupPatchField<volTensorField, tensor>("sigmaExp");


    // Local resistance (defined also for non planestrain cases for consistency)
    scalarField R_P = 1/patch().deltaCoeffs()/(max(2*mu + lambda, VSMALL));

    // Calculate the new gradient R_P*deltaCoeff is equivalent to 1/twoMuLambda
    // for cases with no planeStrain or cases where RP is uniform
    vectorField newGradient(R_P*patch().deltaCoeffs()*(
        (traction_ - pressure_*nCurrent) - (n & sigmaExp)));

    // Force the normal component of the normal gradient (and thus epsilonZZ)
    // to be constant
    if(planeStrain_)
    {
        // Copy only tangent compnent from standard gradient
        newGradient = (newGradient & (I - sqr(n)));

        // Correct RP to consider plane strain approximation
        R_P = scalarField(nFaces, totalArea/gSum(magSf/R_P));

        // Calculate average traction and sigmaExp
        vector tractionEq(1/totalArea*gSum(
            traction_* magSf - pressure_*nCurrent* magSf));

        vector sigmaExpEq(1/totalArea*gSum(n & sigmaExp * magSf));

        newGradient += R_P*patch().deltaCoeffs()*(sqr(n) & (tractionEq - sigmaExpEq));
    } 

    // Resistance network for spring-dashpot system
    if(fixedSpring_)
    {
        // Spring resistance
        scalarField R_S = magSf/max(fixedSpringModulus_/nFacesTotal, SMALL);

        if(planeStrain_)
        {
            R_S = scalarField(nFaces, totalArea/gSum(magSf/R_S));
        }

        // Dashpot resistance
        scalarField R_D = magSf/max(dashpotModulus_/nFacesTotal, SMALL); 

        // Spring-dashpot resistance
        scalarField R_SD = R_S*R_D/(R_S + R_D);

        // Total resistance
        scalarField Rtot = R_P + R_SD;

        // Add dashpot contribution to gradient
        newGradient += R_P*patch().deltaCoeffs()*(*this/R_D);

        // Update weight factor
        w_ = R_SD/Rtot;
    }

    // Relax gradient
    gradient() = (1-relax_)*gradient() + relax_*w_*newGradient;

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void tractionDisplacementFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Lookup the gradient field (D or DD in dependance of the solver used)
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad" + internalField().name()
        );  

    // Face unit normal
    const vectorField n = patch().nf();

    // Cell centers and delta vectors
    vectorField ownCp(patch().boundaryMesh().mesh().C(), patch().faceCells());
    const vectorField delta = patch().Cf() - ownCp;

    // Non-orthogonal correction vectors
    vectorField k = delta - n*(n&delta);

    Field<vector>::operator=
        (
            w_*(patchInternalField()
            + (k & gradField))
            + gradient()/patch().deltaCoeffs()
        );

    fvPatchField<vector>::evaluate();
} 
 

Foam::tmp<Foam::Field<vector>>
Foam::tractionDisplacementFvPatchVectorField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return w_*pTraits<vector>::one;
}


Foam::tmp<Foam::Field<vector>>
Foam::tractionDisplacementFvPatchVectorField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    // Lookup the gradient field (D or DD in dependance of the solver used)
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad" + internalField().name()
        );  

    // Face unit normal
    const vectorField n = patch().nf();

    // Cell centers and delta vectors
    vectorField ownCp(patch().boundaryMesh().mesh().C(), patch().faceCells());
    const vectorField delta = patch().Cf() - ownCp;

    // Non-orthogonal correction vectors
    vectorField k = delta - n*(n&delta);

    return gradient()/this->patch().deltaCoeffs() + w_*(k & gradField);
}


Foam::tmp<Foam::Field<vector>>
Foam::tractionDisplacementFvPatchVectorField::gradientInternalCoeffs() const
{
    return this->patch().deltaCoeffs()*(w_ - 1)*pTraits<vector>::one;
}


Foam::tmp<Foam::Field<vector>>
Foam::tractionDisplacementFvPatchVectorField::gradientBoundaryCoeffs() const
{
    // Lookup the gradient field (D or DD in dependance of the solver used)
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad" + internalField().name()
        );  

    // Face unit normal
    const vectorField n = patch().nf();

    // Cell centers and delta vectors
    vectorField ownCp(patch().boundaryMesh().mesh().C(), patch().faceCells());
    const vectorField delta = patch().Cf() - ownCp;

    // Non-orthogonal correction vectors
    vectorField k = delta - n*(n&delta);

    return gradient() + (w_ - 1)*(k & gradField)*this->patch().deltaCoeffs();
}


void tractionDisplacementFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);

    if (tractionList_.valid())
    {
        os.writeKeyword("tractionList") << nl;
        os << token::BEGIN_BLOCK << nl;
#ifdef OPENFOAMFOUNDATION  
        tractionList_->write(os);
#elif OPENFOAMESI
        tractionList_->writeData(os);
#endif   
        os << token::END_BLOCK << nl;
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
        os.writeKeyword("pressureList") << nl;
        os << token::BEGIN_BLOCK << nl;
#ifdef OPENFOAMFOUNDATION  
        pressureList_->write(os);
#elif OPENFOAMESI
        pressureList_->writeData(os);
#endif
        os << token::END_BLOCK << nl;
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
    writeEntryIfDifferent<scalar>(os, "relax", 1, relax_);
    writeEntryIfDifferent<word>(os, "stressName", "sigma", stressName_);
    writeEntryIfDifferent<bool>(os, "planeStrain", false, planeStrain_);

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

    os.writeEntryIfDifferent<bool>("fixedSpring", false, fixedSpring_);
    if(fixedSpring_)
    {
        os.writeEntry<scalar>("fixedSpringModulus", fixedSpringModulus_);
        os.writeEntry<scalar>("dashpotModulus", dashpotModulus_);
        os.writeEntry("weights", w_);
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
