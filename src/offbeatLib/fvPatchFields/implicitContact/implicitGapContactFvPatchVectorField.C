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
    const fvPatch& patch = this->patch();

    const fvPatchField<scalar>& mu = 
        patch.lookupPatchField<volScalarField, scalar>("mu");

    // Initialize face height with sqrt of area
    scalarField faceHeights(Foam::sqrt(patch.magSf()));

    const pointField& localPoints = patch.patch().localPoints();
    forAll(patch.patch().localFaces(), faceI)
    {
        if(mag(slip_[faceI]) > SMALL)
        {
            const labelList& pointLabelList = patch.patch().localFaces()[faceI];
            scalar maxHeight = -GREAT;
            scalar minHeight = GREAT;

            forAll(pointLabelList, pointI)
            {
                const label& pointIndex = pointLabelList[pointI];
                
                maxHeight = max
                (
                    maxHeight, 
                    localPoints[pointIndex] & (slip_[faceI]/mag(slip_[faceI]))
                );
                
                minHeight = min
                (
                    minHeight, 
                    localPoints[pointIndex] & (slip_[faceI]/mag(slip_[faceI]))
                );
            }

            faceHeights[faceI] = maxHeight - minHeight; 
        }

    }  

    return mu/faceHeights;
}

Foam::tmp<Foam::scalarField> 
Foam::implicitGapContactFvPatchVectorField::gapWidth() const
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

    vectorField Cf = patch.Cf() + totalDispPatch;
                   
    vectorField nbrCf = regionCoupledPatch_.regionCoupledPatch().interpolate
    (
        nbrPatch.Cf() + totalDispNbrPatch
    );

    return( (nbrCf - Cf) & nf);
}    


void Foam::implicitGapContactFvPatchVectorField::correctWeights() 
{
    // Correction factors for non matching surfaces 
    // (e.g. cylinders with a gap)
    const AMIInterpolation& ami =
        regionCoupledPatch_.owner() ?
        regionCoupledPatch_.regionCoupledPatch().AMI() :
        regionCoupledPatch_.nbrPatch().regionCoupledPatch().AMI();

    scalarField correctionFactor = 
        regionCoupledPatch_.owner()?
        ami.srcWeightsSum() : ami.tgtWeightsSum(); 

    // Reference to patch, neighbor patch and neighbor patchField
    const fvPatch& patch = this->patch();
    const fvPatch& nbrPatch = regionCoupledPatch_.neighbFvPatch();
    const regionCoupledOFFBEATFvPatch& nbrRegionCoupledPatch
        = refCast<const regionCoupledOFFBEATFvPatch>(regionCoupledPatch_.nbrPatch());
    const implicitGapContactFvPatchVectorField& nbr = this->neighbour();

    // Relevant patch quantities 
    // (P stands for owner cell/face, N stands for nbr)

    // Patch face areas and normals
    const scalarField& magSfP = patch.magSf();
    const scalarField magSfN = 
        regionCoupledPatch_.regionCoupledPatch().interpolate(nbrPatch.magSf());

    vectorField nP = patch.nf(); 
    vectorField nN = 
        -regionCoupledPatch_.regionCoupledPatch().interpolate(nbrPatch.nf());
    
    // Patch face centers and cell centers
    const vectorField& cfP = patch.Cf(); 
    vectorField cfP_N = 
    nbrRegionCoupledPatch.regionCoupledPatch().interpolate(patch.Cf());
    vectorField cfN = 
        regionCoupledPatch_.regionCoupledPatch().interpolate(nbrPatch.Cf());

    vectorField cpP = 
        vectorField(patch.boundaryMesh().mesh().C(), patch.faceCells());
    vectorField cpN = 
        regionCoupledPatch_.regionCoupledPatch().interpolate
        (
            vectorField(patch.boundaryMesh().mesh().C(), nbrPatch.faceCells())
        );

    // Initial gap correction
    // NOTE: We us nN and not nP so that the projection of cfP on N lands 
    // on N face. If using nP, this projection might land outside N cell
    // for tilted nN/nP normals.
    vectorField initialGapCorrection = ((cfN - cfP) & nN)*nN;
    // vectorField initialGapCorrection = ((cfN - cfP) & nP)*nP;

    // Delta vectors and delta coeffs
    vectorField deltaP = patch.delta();
    vectorField deltaN = (cpN - (cpP + deltaP + initialGapCorrection*(!glued_)));

    scalarField deltaCoeffsP = 1/max((nP & deltaP), SMALL);
    scalarField deltaCoeffsN = 1/max((nN & deltaN), SMALL);

    // Relevant patch fields

    // Calculate two*Mu+Lambda
    const fvPatchField<scalar>& muP = 
        patch.lookupPatchField<volScalarField, scalar>("mu");
    const scalarField muN = 
        regionCoupledPatch_.regionCoupledPatch().interpolate
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>("mu")
        );

    const fvPatchField<scalar>& lambdaP =
        patch.lookupPatchField<volScalarField, scalar>("lambda");
    const scalarField lambdaN = 
        regionCoupledPatch_.regionCoupledPatch().interpolate
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>("lambda")
        );

    const scalarField twoMuLambdaP = (2*muP + lambdaP);
    const scalarField twoMuLambdaN = (2*muN + lambdaN);

    // Patch internal field and nbr internal field
    // NOTE: only taking pnf as a vectorField instead of patchNeighbourField
    // works in multiprocessor 
    vectorField pnif = patchNeighbourField()();
    vectorField pif = patchInternalField()();

    // Nbr patch field
    vectorField pnf = 
        regionCoupledPatch_.regionCoupledPatch().interpolate(this->neighbour());

    // Explicit component of stress tensor
    const tensorField& sigmaExpP =
        patch.lookupPatchField<volTensorField, tensor>("sigmaExp");
    const tensorField sigmaExpN =
        regionCoupledPatch_.regionCoupledPatch().interpolate
        (
            nbrPatch.lookupPatchField<volTensorField, tensor>("sigmaExp")
        );  
    const symmTensorField& sigma =
        patch.lookupPatchField<volSymmTensorField, symmTensor>("sigma");

    // Nbr grad, non orthogonal vectors and correction terms
    const tensorField gradN =
        regionCoupledPatch_.regionCoupledPatch().interpolate
        (
            nbrPatch.lookupPatchField<volTensorField, tensor>("gradD")
        );
    vectorField kCorrN = 
        pos(nP & deltaN)*((nN - deltaN*deltaCoeffsN) & (I - nP*nP));
    vectorField nonOrthN = twoMuLambdaN * (kCorrN & gradN);
    
    // Prepare slave and master DD fields for slip/friction calculation
    vectorField ownDD(patch.size(), Foam::vector::zero);
    vectorField nbrDD(patch.size(), Foam::vector::zero);

    if(db().foundObject<fvMesh>("referenceMesh"))
    {
        const volVectorField& DD =
            db().lookupObject<volVectorField>("DD");

        ownDD = DD.boundaryField()[patch.index()];
        nbrDD = regionCoupledPatch_.regionCoupledPatch().interpolate(
            DD.boundaryField()[nbrPatch.index()]);
    }
    else
    {
        const volVectorField& D =
            db().lookupObject<volVectorField>("D");
        const volVectorField& DOld =
            db().lookupObject<volVectorField>("D_0");

        ownDD = D.boundaryField()[patch.index()] - 
        DOld.boundaryField()[patch.index()];

        nbrDD = regionCoupledPatch_.regionCoupledPatch().interpolate(
            D.boundaryField()[nbrPatch.index()] -
            DOld.boundaryField()[nbrPatch.index()]);   
    }

    slip_ = slip0_ + neg(gapWidth_)*((nbrDD - ownDD) & (I -nP*nP));
    slip_ = neg(gapWidth_)* slip_;

    // Calculate nbr patch stiffness and penalty factors
    const scalarField nbrBoundaryStiff =
        regionCoupledPatch_.regionCoupledPatch().interpolate
        (
            nbr.boundaryStiffness()
        );
    const scalarField nbrBoundaryShearStiff = 
        regionCoupledPatch_.regionCoupledPatch().interpolate
        (
            nbr.boundaryShearStiffness()
        );

    // Calculate local and nbr "resistances"
    const scalarField R_P = 1/twoMuLambdaP/deltaCoeffsP;
    const scalarField R_N = 1/max(twoMuLambdaN, SMALL)/deltaCoeffsN;

    if(method_ == "implicitPenalty")
    {
        penaltyFact_ = neg(gapWidth_)*
            penaltyScaleFact_*min(boundaryStiffness(), nbrBoundaryStiff);
            
        // Calculate normal spring resistances
        scalarField R_S = glued_ ?  magSfP/GREAT : 1/max(penaltyFact_, SMALL);
            
        // Total local resistance (only due to dashpot and owner cell)
        // and coupled resistance
        const scalarField Rtot = R_P + R_N + R_S;

        wC_ = relax_*(1 - R_P/Rtot) + (1 - relax_)*wC_;            
    }
    else if (method_ == "augmentedLagrangian")
    {
        scalarField avgDelta = 0.5*(1/deltaCoeffsP + 1/deltaCoeffsN);
        
        forAll(penaltyFact_, faceI)
        {
            scalar oldValue = wC_[faceI];
            scalar trueWeightC = 1-R_P[faceI]/(R_P[faceI]+R_N[faceI]);
            scalar penTol = max(avgDelta[faceI]*relTolerance_, absTolerance_);
            
            if(glued_)
            {
                wC_[faceI] = (1 - relax_)*oldValue + relax_*trueWeightC;
            }
            else if ((gapWidth_[faceI] > SMALL) or ((sigma[faceI] & nP[faceI]) & nP[faceI]) > 0.0)
            {
                wC_[faceI] = (1 - relax_)*oldValue + relax_*1;
            }
            else if (-gapWidth_[faceI] > penTol)
            {
                wC_[faceI] = (1 - relax_)*oldValue + relax_*trueWeightC;
            }
            else
            {
                scalar alpha = max(-gapWidth_[faceI] / max(penTol, SMALL), 0.0);
                scalar wInterp = (1 - alpha)*1.0 + alpha*trueWeightC;
                wC_[faceI] = (1 - relax_) * oldValue + relax_ * wInterp;
            }
        } 
    }  
    else
    {
        // Print error
    }
    
    // Coupled weights for tangent component
    penaltyFact_t_ = neg(gapWidth_)*
        penaltyScaleFactFriction_*min(boundaryShearStiffness(), nbrBoundaryShearStiff);
    scalarField R_S_t = glued_ ?  magSfP/GREAT : 1/max(penaltyFact_t_, SMALL);
    const scalarField Rtot_t = R_P + R_N + R_S_t;
    wC_t_ = relax_*(1 - R_P/Rtot_t) + (1 - relax_)*wC_t_;

    //- Normal components of the explicit stress tensor
    vectorField Qp = nP & sigmaExpP;
    vectorField Qn = nP & sigmaExpN;

    // Non orthogonal correction to nbr explicit source term.
    // No correction for owner to be consistent to how OpenFOAM handles
    // boundary non-orthogonality
    Qn += nonOrthN;

    vectorField slipCorrection(slip0_.size(), Foam::vector::zero);

    if(db().foundObject<fvMesh>("referenceMesh"))
    {
        slipCorrection = slip0_ & (I -nP*nP); // projection is not necessary
    }
    else
    {

        const vectorField pf0 = 
        patch.lookupPatchField<volVectorField, vector>("D_0");

        const vectorField pnf0 = 
        regionCoupledPatch_.regionCoupledPatch().interpolate
        (
            nbrPatch.lookupPatchField<volVectorField, vector>("D_0")
        );

        slipCorrection = (slip0_ + pf0 - pnf0) & (I -nP*nP);

    }
    

    
    // Estimation of the friction-to-stick ratio
    //------------------------------------------
    // Normal traction estimation
    vectorField traction_n = (1-wC_)/R_P*
        (pnif - pif + initialGapCorrection*(!glued_) + Qp*R_P + Qn*R_N) & (nP*nP);

    // Tangential traction estimation for attached boundaries
    vectorField traction_t = (1-wC_t_)/R_P*
        (pnif - pif + slipCorrection + Qp*R_P + Qn*R_N) & (I - nP*nP);
    
    scalarField magTractionT = 
        min(mag(traction_t), frictionCoeff_*mag(traction_n));

    // Calculation of the friction-to-stick ratio
    // Calculate frictionToStickRatio
    frictionToStickRatio_ = magTractionT/max(mag(traction_t), SMALL);

    // Calculation of the explicit term for the value:
    valueExp_ = -R_P*Qp - R_P*pressure_*nP +
        correctionFactor*
        (
            ((1-wC_)* (initialGapCorrection*(!glued_) + Qp*R_P + Qn*R_N) & (nP*nP)) +
            ((1-wC_t_)*frictionToStickRatio_*
            (slipCorrection + Qp*R_P + Qn*R_N) & (I - nP*nP))
        );

    // Set explicit component of laplacian that is added via manipulateMatrix
    vectorField transverseExp = correctionFactor * deltaCoeffsP *
    (
        (1-wC_)*((pnif-pif)&(nP*nP)) 
        - cmptMultiply((1-wC_)*snGradTransformDiag(),(pnif-pif))
        + (1-wC_t_)*frictionToStickRatio_*((pnif-pif)&(I - nP*nP)) 
        - cmptMultiply(frictionToStickRatio_*(1-wC_t_)*tGradTransformDiag(),(pnif-pif))
    );

    laplacianExp_ = twoMuLambdaP*magSfP*(deltaCoeffsP*valueExp_ + transverseExp);
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
    wC_(patch().weights()),
    wC_t_(patch().weights()),
    penaltyFact_(this->size(), 1.0),
    penaltyFact_t_(this->size(), 1.0),
    valueExp_(this->size(), vector::zero),
    laplacianExp_(this->size(), vector::zero),
    pressure_(p.size(), 0.0),
    interfaceP_(p.size(), 0.0),
    gapWidth_(this->size(), 0),
    gapCalculated_(false),
    frictionToStickRatio_(this->size(), 0),
    slip_(p.size(), Foam::vector::zero),
    slip0_(p.size(), Foam::vector::zero),
    frictionCoeff_(1),
    penaltyScaleFact_(3.5e3),
    penaltyScaleFactFriction_(3.5e3),
    currentTime_(db().time().value()),
    relax_(1),
    glued_(false),
    method_("implicitPenalty"),
    relTolerance_(SMALL),
    absTolerance_(SMALL),
    offset_(0.0),
    coupled_(false)
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
    wC_(patch().weights()),
    wC_t_(patch().weights()),
    penaltyFact_(rcpvf.penaltyFact_),
    penaltyFact_t_(rcpvf.penaltyFact_t_),
    valueExp_(rcpvf.valueExp_),
    laplacianExp_(rcpvf.laplacianExp_),
    pressure_(rcpvf.pressure_),
    interfaceP_(rcpvf.interfaceP_),
    gapWidth_(this->size(), 0),
    gapCalculated_(false),
    frictionToStickRatio_(rcpvf.frictionToStickRatio_),
    slip_(rcpvf.slip_),
    slip0_(rcpvf.slip0_),
    frictionCoeff_(rcpvf.frictionCoeff_),
    penaltyScaleFact_(rcpvf.penaltyScaleFact_),
    penaltyScaleFactFriction_(rcpvf.penaltyScaleFactFriction_),
    currentTime_(rcpvf.currentTime_),
    relax_(rcpvf.relax_),
    glued_(rcpvf.glued_),
    method_(rcpvf.method_),
    relTolerance_(rcpvf.relTolerance_),
    absTolerance_(rcpvf.absTolerance_),
    offset_(rcpvf.offset_),
    coupled_(rcpvf.coupled_)
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
    wC_(patch().size(), 1.0),
    wC_t_(patch().weights()),
    penaltyFact_(this->size(), 0.0),
    penaltyFact_t_(this->size(), 0.0),
    valueExp_(this->size(), vector::zero),
    laplacianExp_(this->size(), vector::zero),
    pressure_(p.size(), 0.0),
    interfaceP_(this->size(), 0),
    gapWidth_(this->size(), 0),
    gapCalculated_(false),
    frictionToStickRatio_(this->size(), 0),
    slip_(p.size(), Foam::vector::zero),
    slip0_(p.size(), Foam::vector::zero),
    frictionCoeff_(dict.lookupOrDefault<scalar>("frictionCoefficient", 0.0)),
    penaltyScaleFact_(0.0),
    penaltyScaleFactFriction_(readScalar(dict.lookup("penaltyFactorFriction"))),
    currentTime_(db().time().value()),
    relax_(dict.lookupOrDefault<scalar>("relax", 1)),
    glued_(dict.lookupOrDefault<bool>("glued", false)),
    method_(dict.lookupOrDefault<word>("method", "implicitPenalty")),
    relTolerance_(dict.lookupOrDefault<scalar>("relativePenetrationTolerance", SMALL)),
    absTolerance_(dict.lookupOrDefault<scalar>("absolutePenetrationTolerance", SMALL)),
    offset_(dict.lookupOrDefault<scalar>("offset", 0)),
    coupled_(false)
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

    if (!(method_ == "implicitPenalty") 
        and !(method_ == "augmentedLagrangian"))
    {
        FatalErrorIn
        (
            "Foam::implicitGapContactFvPatchVectorField::implicitGapContactFvPatchVectorField"
            "("
            "    const fvPatch& p,"
            "    const DimensionedField<vector, volMesh>& iF,"
            "    const dictionary& dict"
            ")"
        )   << "'method' " << method_ << " not recognized"
            << "for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath() << "."
            << "\nAvailable choices are: " 
            << "\n-'implicitPenalty' (default)"
            << "\n-'augmentedLagrangian'" 
            << exit(FatalError);
    }

    if (method_ == "implicitPenalty")
    {
        penaltyScaleFact_ = readScalar(dict.lookup("penaltyFactor"));
    }

    if (dict.found("weights"))
    {
        wC_ = scalarField("weights", dict, p.size());
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
    wC_(patch().weights()),
    wC_t_(patch().weights()),
    penaltyFact_(rcpvf.penaltyFact_),
    penaltyFact_t_(rcpvf.penaltyFact_t_),
    valueExp_(rcpvf.valueExp_),
    laplacianExp_(rcpvf.laplacianExp_),
    pressure_(rcpvf.pressure_),
    interfaceP_(rcpvf.interfaceP_),
    gapWidth_(this->size(), 0),
    gapCalculated_(false),
    frictionToStickRatio_(rcpvf.frictionToStickRatio_),
    slip_(rcpvf.slip_),
    slip0_(rcpvf.slip0_),
    frictionCoeff_(rcpvf.frictionCoeff_),
    penaltyScaleFact_(rcpvf.penaltyScaleFact_),
    penaltyScaleFactFriction_(rcpvf.penaltyScaleFactFriction_),
    currentTime_(rcpvf.currentTime_),
    relax_(rcpvf.relax_),
    glued_(rcpvf.glued_),
    method_(rcpvf.method_),
    relTolerance_(rcpvf.relTolerance_),
    absTolerance_(rcpvf.absTolerance_),
    offset_(rcpvf.offset_),
    coupled_(rcpvf.coupled_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::implicitGapContactFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchVectorField::autoMap(m);

#ifdef OPENFOAMFOUNDATION
    m(wC_, wC_);
    m(wC_t_, wC_t_);
    m(penaltyFact_, penaltyFact_);
    m(penaltyFact_t_, penaltyFact_t_);
    m(valueExp_, valueExp_);
    m(laplacianExp_, laplacianExp_);
    m(pressure_, pressure_);
    m(interfaceP_, interfaceP_);
    m(gapWidth_, gapWidth_);
    m(frictionToStickRatio_, frictionToStickRatio_);
    m(slip_, slip_);
    m(slip0_, slip0_);

#elif OPENFOAMESI    
    wC_.autoMap(m);
    wC_t_.autoMap(m);
    penaltyFact_.autoMap(m);
    penaltyFact_t_.autoMap(m);
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

    wC_.rmap(dmptf.wC_, addr);
    wC_t_.rmap(dmptf.wC_t_, addr);
    penaltyFact_.rmap(dmptf.penaltyFact_, addr);
    penaltyFact_t_.rmap(dmptf.penaltyFact_t_, addr);
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

    // // Lookup the gradient field
    // const tensorField gradField =
    //     patch.lookupPatchField<volTensorField, tensor>
    //     (
    //         "grad" + internalField().name()
    //     );

    // // Face unit normal
    // const vectorField n = patch.nf();  

    // // Cell centers
    // const labelUList& faceCells = patch.faceCells();
    // vectorField ownCp(patch.boundaryMesh().mesh().C(), faceCells);

    // // Delta vectors
    // const vectorField delta = patch.Cf() - ownCp;

    // // Non-orthogonal correction vectors
    // vectorField k = delta - n*(n&delta);

    return patch.deltaCoeffs()*(*this  - patchInternalField());
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

    if(!patch.owner())
    {
        // TODO check that BC parameters are the same on both patches
        // only once and throw error in case
        frictionCoeff_ = nbr.frictionCoeff_;
    }

    // Apply correction
    if (!gapCalculated_)
    {
        // Update the gap width
        gapWidth_ = 
        // 0.5*
        (
            gapWidth() + offset_
          // + patch.regionCoupledPatch().interpolate(nbr.gapWidth())
        );
        // gapWidth_ = 0.5*(gapWidth() + patch.regionCoupledPatch().interpolate(nbr.gapWidth()));  

        // Update the neighbour gap width
        nbr.gapWidth_ = nbrPatch.regionCoupledPatch().interpolate(gapWidth_);     
   
        // pressure_ *= (!glued_)*pos(gapWidth_);

        // Write gapWidth boundary fields
        // TODO move update of gapWidth in mechanics solver
        const_cast<fvPatchScalarField&>
        (
            patch.lookupPatchField<volScalarField, scalar>("gapWidth")
        ) = gapWidth_;   

        const_cast<fvPatchScalarField&>
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>("gapWidth")
        ) = nbr.gapWidth_; 

        // Extract interface pressure (necessray for thermal boundaries)
        correctWeights();
        interfaceP_ = max(-penaltyFact_*gapWidth_, 0.0);
         
        // forAll(interfaceP_, faceI)
        // {
        //     if(gapWidth_[faceI] >= 0)
        //     {  
        //         interfaceP_[faceI]=  0;
        //     }
        // }

        const_cast<fvPatchScalarField&>
        (
            patch.lookupPatchField<volScalarField, scalar>("interfaceP")
        ) = interfaceP_;

        // Interpolate contact pressure on the master side
        nbr.interfaceP_ = nbrPatch.regionCoupledPatch().interpolate(interfaceP_);

        const_cast<fvPatchScalarField&>
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>("interfaceP")
        ) = nbr.interfaceP_;

        gapCalculated_ = true;
        nbr.gapCalculated_ = true;
    }
    else 
    {
        // Make sure gapWidth is always initialized
        if(!gapCalculated_ and db().foundObject<volScalarField>("gapWidth"))
        {
            gapWidth_ = db().lookupObject<volScalarField>("gapWidth").boundaryField()[patch.index()];
        }

        // pressure_ *= (!glued_)*pos(gapWidth_);
        gapCalculated_ = false;
        nbr.gapCalculated_ = false;
        
        correctWeights();
    }

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

    // Normal vectors
    vectorField n(patch.nf()); 
    
    vectorField::operator=
    (
        patchInternalField() + 
        (1.0 - wC_)*correctionFactor*(
            (this->patchNeighbourField() - this->patchInternalField()) & (n*n))
        + (1.0 - wC_t_)*correctionFactor*frictionToStickRatio_*(
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

void Foam::implicitGapContactFvPatchVectorField::
patchNeighbourField(UList<vector>& pnf) const
{
    const Field<vector>& iField = this->primitiveField();
    const labelUList& nbrFaceCells =
        regionCoupledPatch_.neighbFvPatch().faceCells();

    forAll(pnf, facei)
    {
        pnf[facei] = iField[nbrFaceCells[facei]];
    }
}

Foam::tmp<Foam::vectorField>
Foam::implicitGapContactFvPatchVectorField::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    //this function is never called
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

    return pTraits<vector>::one - correctionFactor*
        (
        (1 - wC_)*snGradTransformDiag() + 
        frictionToStickRatio_*(1 - wC_t_)*tGradTransformDiag()
        );
}


Foam::tmp<Foam::vectorField>
Foam::implicitGapContactFvPatchVectorField::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    //this function is never called
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

    return correctionFactor*((1 - wC_)*snGradTransformDiag() + 
        frictionToStickRatio_*(1 - wC_t_)*tGradTransformDiag());
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


    return this->patch().deltaCoeffs()*(
        correctionFactor*(wC_ - 1)*snGradTransformDiag()
        + frictionToStickRatio_*correctionFactor*(wC_t_ - 1)*tGradTransformDiag());
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

    return this->patch().deltaCoeffs()*(
        correctionFactor*(1 - wC_)*snGradTransformDiag() + 
        frictionToStickRatio_*correctionFactor*(1 - wC_t_)*tGradTransformDiag());
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

    if (method_ == "implicitPenalty")
    {
        writeEntry(os, "penaltyFactor", penaltyScaleFact_);
    }

    writeEntry(os, "penaltyFactorFriction", penaltyScaleFactFriction_);

    writeEntryIfDifferent<word>(os, "method", "implicitPenalty", method_);
    
    writeEntryIfDifferent<scalar>(os, "relativePenetrationTolerance", SMALL, relTolerance_);
    writeEntryIfDifferent<scalar>(os, "absolutePenetrationTolerance", SMALL, absTolerance_);

    writeEntry(os, "weights", wC_);
    writeEntry(os, "pressure", pressure_);
    writeEntry(os, "slip", slip_);
    writeEntry(os, "value", *this);

#elif OPENFOAMESI
    os.writeEntryIfDifferent<scalar>("relax", 1, relax_);
    os.writeEntryIfDifferent<bool>("glued", false, glued_);

    os.writeEntry("frictionCoefficient", frictionCoeff_);

    if (method_ == "implicitPenalty")
    {
        os.writeEntry("penaltyFactor", penaltyScaleFact_);
    }
    
    os.writeEntry("penaltyFactorFriction", penaltyScaleFactFriction_);

    os.writeEntryIfDifferent<word>("method", "implicitPenalty", method_);
    os.writeEntryIfDifferent<scalar>("relativePenetrationTolerance", SMALL, relTolerance_);
    os.writeEntryIfDifferent<scalar>("absolutePenetrationTolerance", SMALL, absTolerance_);

    pressure_.writeEntry("pressure", os);
    wC_.writeEntry("weights", os);
    slip_.writeEntry("slip", os);

    this->writeEntry("value", os);
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
