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
#include "implicitGapContactFvPatchVectorField.H"
#include "regionCoupledOFFBEATFvPatch.H"
#include "Time.H"
#include "fvc.H"
#include "AMIInterpolation.H"
#include "PrimitivePatchInterpolation.H"
#include "volPointInterpolation.H"
#include "plane.H"
#include "line.H"
#include "gapGasModel.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Private members  * * * * * * * * * * * * *//
Foam::scalarField 
Foam::implicitGapContactFvPatchVectorField::boundaryStiffness() const
{
    const regionCoupledOFFBEATFvPatch& patch
        = refCast<const regionCoupledOFFBEATFvPatch>(regionCoupledPatch_);

    const fvPatchField<scalar>& threeK = 
        patch.lookupPatchField<volScalarField, scalar>("threeK");
    
    scalarField K = 1.0/3.0*threeK;

    scalarField volumes(patch.size());
            
    forAll(volumes, cellI)
    {
            volumes[cellI] = patch.boundaryMesh().mesh().V()[
                patch.faceCells()[cellI]];
    }

    scalarField faceAreas(patch.size());
            
    forAll(faceAreas, cellI)
    {
            faceAreas[cellI] = patch.magSf()[cellI];
    }  

    return (K)*faceAreas/volumes;

}

Foam::scalarField 
Foam::implicitGapContactFvPatchVectorField::boundaryShearStiffness() const
{
    const regionCoupledOFFBEATFvPatch& patch
        = refCast<const regionCoupledOFFBEATFvPatch>(regionCoupledPatch_);

    const fvPatchField<scalar>& mu = 
        patch.lookupPatchField<volScalarField, scalar>("mu");

    scalarField volumes(patch.size());
            
    forAll(volumes, cellI)
    {
            volumes[cellI] = patch.boundaryMesh().mesh().V()[patch.faceCells()[cellI]];
    }

    scalarField faceAreas(patch.size());
            
    forAll(faceAreas, cellI)
    {
            faceAreas[cellI] = patch.magSf()[cellI];
            // faceAreas[cellI] = sqr(1.0/(2*patch.deltaCoeffs()[cellI]));
    }  

    return (mu)*faceAreas/volumes;
}

Foam::tmp<Foam::scalarField> 
Foam::implicitGapContactFvPatchVectorField::gapWidth() const
{
    const regionCoupledOFFBEATFvPatch& patch
        = refCast<const regionCoupledOFFBEATFvPatch>(regionCoupledPatch_);
    
    const regionCoupledOFFBEATFvPatch& nbrPatch
        = refCast<const regionCoupledOFFBEATFvPatch>(patch.nbrPatch());

    const fvPatchVectorField& totalDispPatch = 
    patch.lookupPatchField<volVectorField, vector>("D"); 
      
    const fvPatchVectorField& totalDispNbrPatch = 
    nbrPatch.lookupPatchField<volVectorField, vector>("D"); 
        
    vectorField nf = patch.Sf() / patch.magSf();

    vectorField Cf = patch.Cf() + totalDispPatch;
                   
    vectorField nbrCf = regionCoupledPatch_.regionCoupledPatch().interpolate(
        nbrPatch.Cf() + totalDispNbrPatch);

    return ((nbrCf - Cf) & nf);
}    


void Foam::implicitGapContactFvPatchVectorField::correctWeights() 
{
    // Correction factors for non matching surfaces (e.g. cylinders with a gap)
    scalarField correctionFactor(this->size(), 1);

    if (regionCoupledPatch_.owner())
    {
        const AMIInterpolation& ami
            = regionCoupledPatch_.regionCoupledPatch().AMI();

        correctionFactor = ami.srcWeightsSum();
    }
    else
    {
        const AMIInterpolation& ami
            = regionCoupledPatch_.nbrPatch().regionCoupledPatch().AMI();

        correctionFactor = ami.tgtWeightsSum();
    }    

    // Reference to patch and neighbor patch
    const fvPatch& patch = this->patch();
    const fvPatch& nbrPatch = regionCoupledPatch_.neighbFvPatch();
    const implicitGapContactFvPatchVectorField& nbr = this->neighbour();


    // Patch quantities 
    // (P stands for owner cell/face, N stands for nbr)

    // Patch face mag
    const scalarField magSfP(patch.magSf());
    const scalarField magSfN(regionCoupledPatch_.regionCoupledPatch().interpolate(
        nbrPatch.magSf()));

    // Owner, neighbour and common normal vectors
    vectorField nP(patch.nf()); 
    vectorField nN(-regionCoupledPatch_.regionCoupledPatch().interpolate(
        nbrPatch.nf())); 
    vectorField n((nP + nN )/mag(nP + nN));

    // Face centers
    vectorField cfP(patch.Cf()); 
    vectorField cfN(regionCoupledPatch_.regionCoupledPatch().interpolate(
        nbrPatch.Cf()));

    // Calculate initial gap to correct Cf position
    vectorField initialGap = (cfN - cfP);

    // Cell centers
    vectorField cpP(patch.boundaryMesh().mesh().C(), patch.faceCells());
    vectorField cpN(regionCoupledPatch_.regionCoupledPatch().interpolate(
        vectorField(patch.boundaryMesh().mesh().C(), nbrPatch.faceCells())));

    // Delta vectors
    vectorField deltaP(cfP - cpP);    
    vectorField deltaN(cpN - cfN);

    // // Alternative deltaN
    // vectorField deltaN(cpN - (cfN));
    // NOTE: nbrDelta is then adjusted to correct for the presence of the gap
    // vectorField deltaN(cpN - cfP);

    // forAll(*this, faceI)
    // {
    //     if(nbrPatch.size() > 0)
    //     {
    //         if(mag(nN[faceI]) > SMALL)
    //         {

    //             plane facePlane(cfN[faceI], n[faceI]);
    //             scalar intersect(facePlane.normalIntersect(
    //                 cfP[faceI], n[faceI]));

    //             vector projectedDist = intersect*n[faceI];

    //             if(intersect >= 0)
    //             {
    //                 deltaN[faceI] -= projectedDist; 
    //             }
    //         }
    //     }
    // }

    // 1/deltaCoeff scalars
    scalarField deltaCoeffsP(1/max((nP & deltaP), SMALL));
    scalarField deltaCoeffsN(1/max((nN & deltaN), SMALL));


    // Patch fields

    // Calculate two*Mu+Lambda
    const fvPatchField<scalar>& muP(
        patch.lookupPatchField<volScalarField, scalar>("mu"));
    const scalarField muN(regionCoupledPatch_.regionCoupledPatch().interpolate(
        nbrPatch.lookupPatchField<volScalarField, scalar>("mu")));

    const fvPatchField<scalar>& lambdaP(
        patch.lookupPatchField<volScalarField, scalar>("lambda"));
    const scalarField lambdaN(regionCoupledPatch_.regionCoupledPatch().interpolate(
        nbrPatch.lookupPatchField<volScalarField, scalar>("lambda")));

    const scalarField twoMuLambdaP = (2*muP + lambdaP);
    const scalarField twoMuLambdaN = (2*muN + lambdaN);

    // Patch neighbour internal field and patch internal field
    // NOTE: only taking pnf as a vectorField instead of patchNeighbourField
    // works in multiprocessor 
    vectorField pnif = patchNeighbourField()();
    vectorField pif = patchInternalField()();

    // Patch oldTime fields, and nbr patch current field
    const vectorField pf0(patch.lookupPatchField<volVectorField, vector>("D_0"));
    const vectorField pnf0(regionCoupledPatch_.regionCoupledPatch().interpolate(
        nbrPatch.lookupPatchField<volVectorField, vector>("D_0")));
    vectorField pnf(regionCoupledPatch_.regionCoupledPatch().interpolate(
        this->neighbour()));

    // Explicit component of stress tensor
    const tensorField sigmaExpP(
        patch.lookupPatchField<volTensorField, tensor>("sigmaExp")); 
    const tensorField sigmaExpN(
        regionCoupledPatch_.regionCoupledPatch().interpolate(
            nbrPatch.lookupPatchField<volTensorField, tensor>("sigmaExp")));  

    // Gradient patchField (necessary for non-orthogonal correction)
    const tensorField gradP(
        patch.lookupPatchField<volTensorField, tensor>("gradD")); 
    const tensorField gradN(
        regionCoupledPatch_.regionCoupledPatch().interpolate(
            nbrPatch.lookupPatchField<volTensorField, tensor>("gradD")));

    //- Normal components of the explicit stress tensor
    vectorField Qp(nP & sigmaExpP);
    vectorField Qn(nN & sigmaExpN);


    // Non orthogonal correction

    // Non orthogonal vectors
    vectorField kCorrP(pos(nP & deltaP)*(nP - (deltaP*deltaCoeffsP)));
    vectorField kCorrN(pos(nN & deltaN)*(nN - (deltaN*deltaCoeffsN)));
    
    // Non orthogonal correction terms
    vectorField nonOrthP = twoMuLambdaP * (kCorrP & gradP);
    vectorField nonOrthN = twoMuLambdaN * (kCorrN & gradN);

    // Add non orthogonal correction to explicit source term.
    // Add only normal component, as for the tangent component this correction
    // is included in the explicit gradient (i.e. like a fixedGradient BC)
    Qp += nonOrthP;
    Qn += nonOrthN;

    // Increase slip (only when penetration is > 0)
    // When gapWidth is zero total slip is equal to zero
    slip_ = slip0_ + neg(gapWidth_)*(
        ((pnf - pnf0) - (*this - pf0)) & (I -nP*nP));
    slip_ = neg(gapWidth_)* slip_;

    // Calculate "resistances" network

    // Interpolate back to this patch, so to avoid scalarFields of different lenght
    const scalarField nbrBoundaryStiff(
        regionCoupledPatch_.regionCoupledPatch().interpolate(nbr.boundaryStiffness()));
    const scalarField nbrBoundaryShearStiff(
        regionCoupledPatch_.regionCoupledPatch().interpolate(nbr.boundaryShearStiffness()));
    
    scalarField avgDelta(0.5*(1/deltaCoeffsP/10 + 1/deltaCoeffsN/10));

     scalarField newPenaltyFact = pow(min(1/(avgDelta)*
        max(-gapWidth_, SMALL), 1.0)/0.03, 2.0)/4.0*0.03*penaltyScaleFact_*min(
            boundaryStiffness(), nbrBoundaryStiff);

    scalarField newPenaltyFact_dashpot = min(1/(avgDelta)*
        max(-gapWidth_, SMALL), 1.0)*penaltyScaleFactDashpot_*min(
            boundaryStiffness(), nbrBoundaryStiff);

    scalarField newPenaltyFact_t = min(1/(avgDelta)*
        max(-gapWidth_, SMALL), 1.0)*penaltyScaleFactFriction_*min(
            boundaryShearStiffness(), nbrBoundaryShearStiff);

    // Spring resistance
    penaltyFact_ = (
        relax_*newPenaltyFact + (1 - relax_)*penaltyFact_);
    scalarField R_S = glued_ ?  magSfP/GREAT : 1/max(penaltyFact_, SMALL);

    // Spring resistance for tangent component
    penaltyFact_t_ = (
        relax_*newPenaltyFact_t + (1 - relax_)*penaltyFact_t_);
    scalarField R_S_t = glued_ ?  magSfP/GREAT : 1/max(penaltyFact_t_, SMALL);

    // Dashpot resistance
    penaltyFactDashpot_ = (
        relax_*newPenaltyFact_dashpot + (1 - relax_)*penaltyFactDashpot_);
    scalarField R_D = glued_ ?  magSfP/SMALL : 1/max(penaltyFactDashpot_, SMALL);

    // Calculate normal "resistances"
    const scalarField R_P = 1/twoMuLambdaP/deltaCoeffsP;
    const scalarField R_N = 1/max(twoMuLambdaN, SMALL)/deltaCoeffsN;

    // Total local resistance (only due to dashpot and owner cell)
    const scalarField RtotP = R_P + R_D;

    // Total coupled resistance (due to spring, nbr and owner cell)
    const scalarField RtotC = R_P + R_N + R_S;

    // Local weights
    wP_ = R_D/RtotP;

    // Coupled weights
    wC_ = 1 - R_P/RtotC;

    // Coupled weights for tangent component
    wC_t_ = 1 - R_P/(R_P + R_N + R_S_t);
      

    // Initial gap correction

    if(!glued_)
    {
        vectorField initialGapCorrection(initialGap/R_N);
        Qn += (initialGapCorrection & nP)*nP;
    }   
      

    // Slip correction

    if(!glued_)
    {
        vectorField slipCorrection((slip0_ + pf0 - pnf0)/R_N);
        Qn += (slipCorrection & (I -nP*nP));
    }    

    Qn *= correctionFactor;

    // Set explicit components of gradient, value and laplacian vectors

    // Set gradient explicit components due to contact and pressure-dashpot
    vectorField contactGradExp_n(deltaCoeffsP*(1 - wC_)*(
            (Qp*R_P + Qn*R_N + 0*pressure_*(R_P + R_N)*nP) & (nP*nP)));
    vectorField contactGradExp_t(deltaCoeffsP*(1 - wC_t_)*(
            (Qp*R_P + Qn*R_N + 0*pressure_*(R_P + R_N)*nP) & (I - nP*nP)));
    vectorField pressreDashpotGradExp(
        (-pressure_*nP - Qp + *this/R_D)/twoMuLambdaP);

    // Full traction (trial/explicit for determining maximum tangent force)
    vectorField traction_n( twoMuLambdaP*(
        deltaCoeffsP*(1 - wC_)*((pnif - pif) & (nP*nP))+ contactGradExp_n));
    vectorField traction_t( twoMuLambdaP*(
        deltaCoeffsP*(1 - wC_t_)*((pnif - pif) & (I -nP*nP))+ contactGradExp_t));

    // Apply friction limiter for tangential traction.
    scalarField magTractionT = min(
        mag(traction_t), frictionCoeff_*mag(traction_n));

    // Calculate frictionToStickRatio
    frictionToStickRatio_ = magTractionT/max(mag(traction_t), SMALL);

    // Adjust explicit gradient adding normal and tangent components
    gradExp_ = wP_*(pressreDashpotGradExp + contactGradExp_n
        + frictionToStickRatio_*contactGradExp_t);

    // Set explicit component of value
    valueExp_ = gradExp_/deltaCoeffsP;

    // Set explicit component of laplacian that is added via manipulateMatrix
    laplacianExp_ = (gradExp_ + (kCorrP & gradP))*twoMuLambdaP*magSfP;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::implicitGapContactFvPatchVectorField::
implicitGapContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    coupledFvPatchField<vector>(p, iF),
    regionCoupledPatch_(refCast<const regionCoupledBaseOFFBEATFvPatch>(p)),  
    wP_(patch().weights()),
    wC_(patch().weights()),
    wC_t_(patch().weights()),
    penaltyFact_(this->size(), 1.0),
    penaltyFact_t_(this->size(), 1.0),
    penaltyFactDashpot_(this->size(), 1.0),    
    gradExp_(this->size(), vector::zero),
    valueExp_(this->size(), vector::zero),
    laplacianExp_(this->size(), vector::zero),
    pressure_(p.size(), 0.0),
    interfaceP_(p.size(), 0.0),
    gapWidth_(this->size(), 0),
    frictionToStickRatio_(this->size(), 0),
    slip_(p.size(), Foam::vector::zero),
    slip0_(p.size(), Foam::vector::zero),
    frictionCoeff_(1),
    penaltyScaleFact_(3.5e3),
    penaltyScaleFactFriction_(3.5e3),
    penaltyScaleFactDashpot_(3.5e3),
    currentTime_(db().time().value()),
    relax_(1),
    glued_(false)
{}


Foam::implicitGapContactFvPatchVectorField::
implicitGapContactFvPatchVectorField
(
    const implicitGapContactFvPatchVectorField& rcpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<vector>(rcpvf, p, iF, mapper),
    regionCoupledPatch_(refCast<const regionCoupledBaseOFFBEATFvPatch>(p)),
    wP_(patch().weights()),
    wC_(patch().weights()),
    wC_t_(patch().weights()),
    penaltyFact_(rcpvf.penaltyFact_),
    penaltyFact_t_(rcpvf.penaltyFact_t_),
    penaltyFactDashpot_(rcpvf.penaltyFact_t_),
    gradExp_(rcpvf.gradExp_),
    valueExp_(rcpvf.valueExp_),
    laplacianExp_(rcpvf.laplacianExp_),
    pressure_(rcpvf.pressure_),
    interfaceP_(rcpvf.interfaceP_),
    gapWidth_(this->size(), 0),
    frictionToStickRatio_(rcpvf.frictionToStickRatio_),
    slip_(rcpvf.slip_),
    slip0_(rcpvf.slip0_),
    frictionCoeff_(rcpvf.frictionCoeff_),
    penaltyScaleFact_(rcpvf.penaltyScaleFact_),
    penaltyScaleFactFriction_(rcpvf.penaltyScaleFactFriction_),
    penaltyScaleFactDashpot_(rcpvf.penaltyScaleFactDashpot_),
    currentTime_(rcpvf.currentTime_),
    relax_(rcpvf.relax_),
    glued_(rcpvf.glued_)
{}


Foam::implicitGapContactFvPatchVectorField::
implicitGapContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<vector>(p, iF, dict),
    regionCoupledPatch_(refCast<const regionCoupledBaseOFFBEATFvPatch>(p)),
    wP_(patch().weights()),
    wC_(patch().weights()),
    wC_t_(patch().weights()),
    penaltyFact_(this->size(), 0.0),
    penaltyFact_t_(this->size(), 0.0),
    penaltyFactDashpot_(this->size(), 0.0),
    gradExp_(this->size(), vector::zero),
    valueExp_(this->size(), vector::zero),
    laplacianExp_(this->size(), vector::zero),
    pressure_(p.size(), 0.0),
    interfaceP_(this->size(), 0),
    gapWidth_(this->size(), 0),
    frictionToStickRatio_(this->size(), 0),
    slip_(p.size(), Foam::vector::zero),
    slip0_(p.size(), Foam::vector::zero),
    frictionCoeff_(dict.lookupOrDefault<scalar>("frictionCoefficient", 0.0)),
    penaltyScaleFact_(readScalar(dict.lookup("penaltyFactor"))),
    penaltyScaleFactFriction_(readScalar(dict.lookup("penaltyFactorFriction"))),
    penaltyScaleFactDashpot_(readScalar(dict.lookup("penaltyFactorDashpot"))),
    currentTime_(db().time().value()),
    relax_(dict.lookupOrDefault<scalar>("relax", 1)),
    glued_(dict.lookupOrDefault<bool>("glued", false))
{    
    if (!isA<regionCoupledBaseOFFBEAT>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "Foam::implicitGapContactFvPatchVectorField::implicitGapContactFvPatchVectorField"
            "("
            "    const fvPatch& p,"
            "    const DimensionedField<vector, volMesh>& iF,"
            "    const dictionary& dict"
            ")"
        )   << "' not type '" << regionCoupledBaseOFFBEAT::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    if (dict.found("pressure"))
    {
        pressure_ = scalarField("pressure", dict, patch().size());
    }     

    if (dict.found("slip"))
    {
        slip_ = vectorField("slip", dict, p.size());
        slip0_ = slip_;
    }
}


Foam::implicitGapContactFvPatchVectorField::
implicitGapContactFvPatchVectorField
(
    const implicitGapContactFvPatchVectorField& rcpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    coupledFvPatchField<vector>(rcpvf, iF),
    regionCoupledPatch_(rcpvf.regionCoupledPatch_),
    wP_(patch().weights()),
    wC_(patch().weights()),
    wC_t_(patch().weights()),
    penaltyFact_(rcpvf.penaltyFact_),
    penaltyFact_t_(rcpvf.penaltyFact_t_),
    penaltyFactDashpot_(rcpvf.penaltyFact_t_),
    gradExp_(rcpvf.gradExp_),
    valueExp_(rcpvf.valueExp_),
    laplacianExp_(rcpvf.laplacianExp_),
    pressure_(rcpvf.pressure_),
    interfaceP_(rcpvf.interfaceP_),
    gapWidth_(this->size(), 0),
    frictionToStickRatio_(rcpvf.frictionToStickRatio_),
    slip_(rcpvf.slip_),
    slip0_(rcpvf.slip0_),
    frictionCoeff_(rcpvf.frictionCoeff_),
    penaltyScaleFact_(rcpvf.penaltyScaleFact_),
    penaltyScaleFactFriction_(rcpvf.penaltyScaleFactFriction_),
    penaltyScaleFactDashpot_(rcpvf.penaltyScaleFactDashpot_),
    currentTime_(rcpvf.currentTime_),
    relax_(rcpvf.relax_),
    glued_(rcpvf.glued_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::implicitGapContactFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchVectorField::autoMap(m);

#ifdef OPENFOAMFOUNDATION
    m(wP_, wP_);
    m(wC_, wC_);
    m(wC_t_, wC_t_);
    m(penaltyFact_, penaltyFact_);
    m(penaltyFact_t_, penaltyFact_t_);
    m(penaltyFactDashpot_, penaltyFactDashpot_);
    m(gradExp_, gradExp_);
    m(valueExp_, valueExp_);
    m(laplacianExp_, laplacianExp_);
    m(pressure_, pressure_);
    m(interfaceP_, interfaceP_);
    m(gapWidth_, gapWidth_);
    m(frictionToStickRatio_, frictionToStickRatio_);
    m(slip_, slip_);
    m(slip0_, slip0_);

#elif OPENFOAMESI    
    wP_.autoMap(m);
    wC_.autoMap(m);
    wC_t_.autoMap(m);
    penaltyFact_.autoMap(m);
    penaltyFact_t_.autoMap(m);
    penaltyFactDashpot_.autoMap(m);
    gradExp_.autoMap(m);
    valueExp_.autoMap(m);
    laplacianExp_.autoMap(m);
    pressure_.autoMap(m);
    interfaceP_.autoMap(m);
    gapWidth_.autoMap(m);
    frictionToStickRatio_.autoMap(m);
    slip_.autoMap(m);
    slip0_.autoMap(m);    
#endif    
}


void Foam::implicitGapContactFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fvPatchVectorField::rmap(ptf, addr);

    const implicitGapContactFvPatchVectorField& dmptf =
        refCast<const implicitGapContactFvPatchVectorField>(ptf);

    wP_.rmap(dmptf.wP_, addr);
    wC_.rmap(dmptf.wC_, addr);
    wC_t_.rmap(dmptf.wC_t_, addr);
    penaltyFact_.rmap(dmptf.penaltyFact_, addr);
    penaltyFact_t_.rmap(dmptf.penaltyFact_t_, addr);
    penaltyFactDashpot_.rmap(dmptf.penaltyFactDashpot_, addr);
    gradExp_.rmap(dmptf.gradExp_, addr);
    valueExp_.rmap(dmptf.valueExp_, addr);
    laplacianExp_.rmap(dmptf.laplacianExp_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
    interfaceP_.rmap(dmptf.interfaceP_, addr);
    gapWidth_.rmap(dmptf.gapWidth_, addr);
    frictionToStickRatio_.rmap(dmptf.frictionToStickRatio_, addr);
}


Foam::tmp<Foam::vectorField> Foam::implicitGapContactFvPatchVectorField::
snGrad() const
{
    const fvPatch& patch = this->patch();

    // Lookup the gradient field
    const tensorField gradField =
        patch.lookupPatchField<volTensorField, tensor>
        (
            "grad" + internalField().name()
        );

    // Face unit normal
    const vectorField n = patch.nf();  

    // Cell centers
    const labelUList& faceCells = patch.faceCells();
    vectorField ownCp(patch.boundaryMesh().mesh().C(), faceCells);

    // Delta vectors
    const vectorField delta = patch.Cf() - ownCp;

    // Non-orthogonal correction vectors
    vectorField k = delta - n*(n&delta);

    return patch.deltaCoeffs()*(*this - (k & gradField) - patchInternalField());
}


Foam::tmp<Foam::vectorField> Foam::implicitGapContactFvPatchVectorField::
snGrad(const scalarField&) const
{
    return snGrad();
}


void Foam::implicitGapContactFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Update gap pressure
    pressure_ = db().lookupObject<gapGasModel>("gapGas").p();

    const regionCoupledOFFBEATFvPatch& patch
        = refCast<const regionCoupledOFFBEATFvPatch>(regionCoupledPatch_);
    
    const regionCoupledOFFBEATFvPatch& nbrPatch
        = refCast<const regionCoupledOFFBEATFvPatch>(patch.nbrPatch());

    const implicitGapContactFvPatchVectorField& nbr = this->neighbour();

    if(currentTime_ < db().time().value())
    {
        slip0_ = slip_;
        currentTime_ = db().time().value();
    }

    //- Apply correction
    if (!regionCoupledPatch_.owner())
    {
        frictionCoeff_ = nbr.frictionCoeff_;

        const AMIInterpolation& ami(nbrPatch.regionCoupledPatch().AMI());

        // Update the gap width
        gapWidth_ = 0.5*(gapWidth() + 
            regionCoupledPatch_.regionCoupledPatch().interpolate(nbr.gapWidth()));    

        // Debugging
        const_cast<fvPatchScalarField&>
        (
            patch.lookupPatchField<volScalarField, scalar>("gapWidth")
        ) = gapWidth_;    

        // The explicit terms are lumped together in Q
        const symmTensorField stress =
            patch.lookupPatchField<volSymmTensorField, symmTensor>("sigma"); 

        vectorField n(patch.nf() );

        interfaceP_ = -( ((stress&n)&n) );
         
        forAll(interfaceP_, faceI)
        {
            if(gapWidth_[faceI] >= 0)
            {  
                interfaceP_[faceI]=  0;
            }
            else
            {
                interfaceP_[faceI]= max(interfaceP_[faceI] - pressure_[faceI], 0.0); 
            }
        }

        const_cast<fvPatchScalarField&>
        (
            patch.lookupPatchField<volScalarField, scalar>("interfaceP")
        ) = interfaceP_;
    }
    else 
    {
        const AMIInterpolation& ami(regionCoupledPatch_.regionCoupledPatch().AMI());

        // Update the gap width
        gapWidth_ = 0.5*(gapWidth() + 
            regionCoupledPatch_.regionCoupledPatch().interpolate(nbr.gapWidth()));   

        // Debugging
        const_cast<fvPatchScalarField&>
        (
            patch.lookupPatchField<volScalarField, scalar>("gapWidth")
        ) = gapWidth_; 

        // The explicit terms are lumped together in Q
        const symmTensorField stress =
            patch.lookupPatchField<volSymmTensorField, symmTensor>("sigma"); 

        vectorField ownN(patch.nf()); 
        vectorField nbrN(-regionCoupledPatch_.regionCoupledPatch().interpolate(
            regionCoupledPatch_.neighbFvPatch().nf()));  

        vectorField n(patch.nf() );

        interfaceP_ = -( ((stress&n)&n) );
         
        forAll(interfaceP_, faceI)
        {
            if(gapWidth_[faceI] >= 0)
            {  
                interfaceP_[faceI]=  0;
            }
            else
            {
                interfaceP_[faceI]= max(interfaceP_[faceI] - pressure_[faceI], 0.0); 
            }
        }

        const_cast<fvPatchScalarField&>
        (
            patch.lookupPatchField<volScalarField, scalar>("interfaceP")
        ) = interfaceP_;
    }

    correctWeights();

    fvPatchVectorField::updateCoeffs();
}


void Foam::implicitGapContactFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!updated())
    {
        updateCoeffs();
    } 

    // Reference to patch and neighbor patch
    const fvPatch& patch = this->patch();
    const fvPatch& nbrPatch = regionCoupledPatch_.neighbFvPatch();

    // Patch quantities 
    // (P stands for owner cell/face, N stands for nbr)

    // Owner, neighbour and common normal vectors
    vectorField n(patch.nf()); 
    
    vectorField::operator=
    (
        wP_*patchInternalField()
        + wP_*(1.0 - wC_)*(
            (this->patchNeighbourField() - this->patchInternalField()) & (n*n))
        + wP_*(1.0 - wC_t_)*frictionToStickRatio_*(
            (this->patchNeighbourField() - this->patchInternalField()) & (I - n*n))
        + valueExp_
    ); 

    fvPatchVectorField::evaluate();
}


Foam::tmp<Foam::Field<Foam::vector> >
Foam::implicitGapContactFvPatchVectorField::
patchNeighbourField() const
{
    const fvPatch& nbrPatch = regionCoupledPatch_.neighbFvPatch();

    const labelUList& nbrFaceCells = nbrPatch.faceCells();
    
    const vectorField npf(internalField(), nbrFaceCells);

    return regionCoupledPatch_.regionCoupledPatch().interpolate(npf);
}

Foam::tmp<Foam::vectorField>
Foam::implicitGapContactFvPatchVectorField::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    Foam::vector ones(pTraits<vector>::one);

    tmp<vectorField> vICoeffs
    (
        new vectorField(this->size(), ones)
    ); 
    vectorField& vICoeffsRef = vICoeffs.ref();

    vICoeffsRef *= wP_;

    return vICoeffs;
}


Foam::tmp<Foam::vectorField>
Foam::implicitGapContactFvPatchVectorField::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    Foam::vector ones(pTraits<vector>::one);

    tmp<vectorField> vBCoeffs
    (
        new vectorField(this->size(), ones)
    ); 
    vectorField& vBCoeffsRef = vBCoeffs.ref();

    vBCoeffsRef *= (1 - wP_);

    return vBCoeffs;
}


Foam::tmp<Foam::vectorField>
Foam::implicitGapContactFvPatchVectorField::gradientInternalCoeffs
(
    const scalarField&
) const
{
    return gradientInternalCoeffs();
}


Foam::tmp<Foam::vectorField>
Foam::implicitGapContactFvPatchVectorField::gradientInternalCoeffs() const
{ 
    return this->patch().deltaCoeffs()*(
        (wP_ - 1)*pTraits<vector>::one + wP_*(wC_ - 1)*snGradTransformDiag()
        + frictionToStickRatio_*wP_*(wC_t_ - 1)*tGradTransformDiag());
}


Foam::tmp<Foam::vectorField>
Foam::implicitGapContactFvPatchVectorField::gradientBoundaryCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    return gradientBoundaryCoeffs();
}


Foam::tmp<Foam::vectorField>
Foam::implicitGapContactFvPatchVectorField::gradientBoundaryCoeffs() const
{  
    return this->patch().deltaCoeffs()*(
        wP_*(1 - wC_)*snGradTransformDiag() + 
        frictionToStickRatio_*wP_*(1 - wC_t_)*tGradTransformDiag());
}


#ifdef OPENFOAMFOUNDATION
void Foam::implicitGapContactFvPatchVectorField::updateInterfaceMatrix
(
    Field<scalar>& result,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes comms
) const
{
    const labelUList& faceCells = regionCoupledPatch_.faceCells();
    const labelUList& nbrFaceCells
        = regionCoupledPatch_.neighbFvPatch().faceCells();
    
    const scalarField pnf
        = regionCoupledPatch_.regionCoupledPatch().interpolate
        (
            scalarField(psiInternal, nbrFaceCells)
        );

    // Multiply the field by coefficients and add into the result
    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
    }
}


void Foam::implicitGapContactFvPatchVectorField::updateInterfaceMatrix
(
    Field<vector>& result,
    const Field<vector>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes comms
) const
{
    const labelUList& faceCells = regionCoupledPatch_.faceCells();
    const labelUList& nbrFaceCells
        = regionCoupledPatch_.neighbFvPatch().faceCells();
    
    const vectorField pnf
        = regionCoupledPatch_.regionCoupledPatch().interpolate
        (
            vectorField(psiInternal, nbrFaceCells)
        );

    // Multiply the field by coefficients and add into the result
    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
    }
}

#elif OPENFOAMESI
void Foam::implicitGapContactFvPatchVectorField::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    const labelUList& faceCells = regionCoupledPatch_.faceCells();
    const labelUList& nbrFaceCells
        = regionCoupledPatch_.neighbFvPatch().faceCells();
 
    solveScalarField pnf
        = regionCoupledPatch_.regionCoupledPatch().interpolate
        (
            scalarField(psiInternal, nbrFaceCells)
        ); 
 
    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}


void Foam::implicitGapContactFvPatchVectorField::updateInterfaceMatrix
(
    Field<vector>& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<vector>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{
    const labelUList& faceCells = regionCoupledPatch_.faceCells();
    const labelUList& nbrFaceCells
        = regionCoupledPatch_.neighbFvPatch().faceCells();

    const vectorField pnf
        = regionCoupledPatch_.regionCoupledPatch().interpolate
        (
            vectorField(psiInternal, nbrFaceCells)
        );

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}
#endif


void Foam::implicitGapContactFvPatchVectorField::manipulateMatrix(fvMatrix<vector>& matrix)
{
    const labelUList& faceCells = regionCoupledPatch_.faceCells();

    //- Create field to be added to the matrix source
    //- It is the explicit component of the Laplacian
    //- equal to the explicit component of the gradient
    // - multiplied by K and magnitude of face
    vectorField source(matrix.source().size(), Foam::vector(0, 0, 0));

    forAll(faceCells, elemI)
    {
        source[faceCells[elemI]] = laplacianExp_[elemI];
    }

    matrix.source() += source;
} 


Foam::tmp<Foam::Field<Foam::vector>>
Foam::implicitGapContactFvPatchVectorField::snGradTransformDiag() const
{
    symmTensorField valueFraction(sqr(patch().nf()));
    vectorField diag(valueFraction.size());

    diag.replace
    (
        vector::X,
        sqrt(mag(valueFraction.component(symmTensor::XX)))
    );
    diag.replace
    (
        vector::Y,
        sqrt(mag(valueFraction.component(symmTensor::YY)))
    );
    diag.replace
    (
        vector::Z,
        sqrt(mag(valueFraction.component(symmTensor::ZZ)))
    );

 return transformFieldMask<Foam::vector>(
    pow<vector, pTraits<Foam::vector>::rank>(diag));
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::implicitGapContactFvPatchVectorField::tGradTransformDiag() const
{
    symmTensorField valueFraction(I - sqr(patch().nf()));
    vectorField diag(valueFraction.size());

    diag.replace
    (
        vector::X,
        sqrt(mag(valueFraction.component(symmTensor::XX)))
    );
    diag.replace
    (
        vector::Y,
        sqrt(mag(valueFraction.component(symmTensor::YY)))
    );
    diag.replace
    (
        vector::Z,
        sqrt(mag(valueFraction.component(symmTensor::ZZ)))
    );

 return transformFieldMask<Foam::vector>(
    pow<vector, pTraits<Foam::vector>::rank>(diag));
}


void Foam::implicitGapContactFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);

#ifdef OPENFOAMFOUNDATION  
    writeEntryIfDifferent<scalar>(os, "relax", 1, relax_);
    writeEntryIfDifferent<bool>(os, "glued", false, glued_);

    writeEntry(os, "frictionCoefficient", frictionCoeff_);
    writeEntry(os, "penaltyFactor", penaltyScaleFact_);
    writeEntry(os, "penaltyFactorFriction", penaltyScaleFactFriction_);
    writeEntry(os, "penaltyFactorDashpot", penaltyScaleFactDashpot_);

    writeEntry(os, "pressure", pressure_);
    writeEntry(os, "slip", slip_);
    writeEntry(os, "value", *this);

#elif OPENFOAMESI
    os.writeEntryIfDifferent<scalar>("relax", 1, relax_);
    os.writeEntryIfDifferent<bool>("glued", false, glued_);

    os.writeEntry("frictionCoefficient", frictionCoeff_);
    os.writeEntry("penaltyFactor", penaltyScaleFact_);
    os.writeEntry("penaltyFactorFriction", penaltyScaleFactFriction_);
    os.writeEntry("penaltyFactorDashpot", penaltyScaleFactDashpot_);

    os.writeEntry<scalarField>("pressure", pressure_);
    os.writeEntry("slip", slip_);

    writeEntry("value", os);
#endif    
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        implicitGapContactFvPatchVectorField
    );
};


// ************************************************************************* //
