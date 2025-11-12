/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2506                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2025 OpenCFD Ltd.          |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "junctionJumpAMIFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "primitiveMeshTools.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::junctionJumpAMIFvPatchVectorField::junctionJumpAMIFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    jumpDiscontinuousCyclicAMIFvPatchField<vector>(p, iF),
    jump_(this->size(), Zero),
    residual_(Zero),
    residualPrev_(Zero),
    underRelaxation_(1),
    aitkenRelaxation_(0.01),
    nIter_(0),
    useAitken_(false),
    scalingFactors_(this->size(), Zero),
    overlapPatchNames_(this->size(), ""),
    startTime_(0)
{}


Foam::junctionJumpAMIFvPatchVectorField::junctionJumpAMIFvPatchVectorField
(
    const junctionJumpAMIFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    jumpDiscontinuousCyclicAMIFvPatchField<vector>(ptf, p, iF, mapper),
    jump_(ptf.jump_, mapper),
    residual_(ptf.residual_),
    residualPrev_(ptf.residualPrev_),
    underRelaxation_(ptf.underRelaxation_),
    aitkenRelaxation_(ptf.aitkenRelaxation_),
    nIter_(ptf.nIter_),
    useAitken_(ptf.useAitken_),
    scalingFactors_(ptf.scalingFactors_),
    overlapPatchNames_(ptf.overlapPatchNames_),
    startTime_(ptf.startTime_)
{}


Foam::junctionJumpAMIFvPatchVectorField::junctionJumpAMIFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    jumpDiscontinuousCyclicAMIFvPatchField<vector>(p, iF),
    jump_(p.size(), Zero),
    residual_(Zero),
    residualPrev_(Zero),
    underRelaxation_(dict.getOrDefault<scalar>("underRelaxation", 1)),
    aitkenRelaxation_(dict.getOrDefault<scalar>("aitkenRelaxation", 0.01)),
    nIter_(0),
    useAitken_(false),
    scalingFactors_(this->size(), Zero),
    overlapPatchNames_(this->size(), ""),
    startTime_(dict.getOrDefault<scalar>("startTime", 0))
{

    if (this->discontinuousCyclicAMIPatch().owner())
    {
        jump_.assign("jump", dict, p.size(), IOobjectOption::LAZY_READ);
    }


    if(!this->discontinuousCyclicAMIPatch().owner())
    {
        overlapPatchNames_ = dict.get<wordList>("connectedPatches");
    }

    if(dict.found("aitkenAcceleration"))
        useAitken_ = dict.get<bool>("aitkenAcceleration");


    // Modify areas

    if(this->discontinuousCyclicAMIPatch().owner())
    {
        scalar nBranch = overlapPatchNames_.size()+1;

        // Rescale polypatch areas

        discontinuousCyclicAMIPolyPatch& discPP = 
            const_cast<discontinuousCyclicAMIPolyPatch&>
            (
                this->discontinuousCyclicAMIPatch().cyclicAMIPatch()
            );

        vectorField::subField Sf = discPP.faceAreas();
        forAll(Sf, facei)
        {
            Sf[facei] = Sf[facei] / nBranch;
        }

        discPP.areaFraction(1/nBranch);

        const polyMesh& mesh = this->patch().boundaryMesh().mesh();


        primitiveMeshTools::updateCellCentresAndVols
        (
            mesh,
            mesh.faceCentres(),
            mesh.faceAreas(),                      
            uniqueSort(discPP.faceCells()), 
            mesh.cells(),
            const_cast<vectorField&>(mesh.cellCentres()),
            const_cast<scalarField&>(mesh.cellVolumes())
        );

        // // Also modify fvPatch areas

        const discontinuousCyclicAMIFvPatch& discFVP = 
        (
            this->discontinuousCyclicAMIPatch()
        );

        const_cast<vectorField&>(discFVP.Sf()) = discFVP.patch().faceAreas();
        const_cast<vectorField&>(discFVP.Cf()) = discFVP.patch().faceCentres();
        const_cast<scalarField&>(discFVP.magSf()) = mag(discFVP.patch().faceAreas());


        this->patch().boundaryMesh().mesh().V().write();
    }
        
}



Foam::junctionJumpAMIFvPatchVectorField::junctionJumpAMIFvPatchVectorField
(
    const junctionJumpAMIFvPatchVectorField& ptf
)
:
    jumpDiscontinuousCyclicAMIFvPatchField<vector>(ptf),
    jump_(ptf.jump_),
    residual_(ptf.residual_),
    residualPrev_(ptf.residualPrev_),
    underRelaxation_(ptf.underRelaxation_),
    aitkenRelaxation_(ptf.aitkenRelaxation_),
    nIter_(ptf.nIter_),
    useAitken_(ptf.useAitken_),
    scalingFactors_(ptf.scalingFactors_),
    overlapPatchNames_(ptf.overlapPatchNames_),
    startTime_(ptf.startTime_)
{}


Foam::junctionJumpAMIFvPatchVectorField::junctionJumpAMIFvPatchVectorField
(
    const junctionJumpAMIFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    jumpDiscontinuousCyclicAMIFvPatchField<vector>(ptf, iF),
    jump_(ptf.jump_),
    residual_(ptf.residual_),
    residualPrev_(ptf.residualPrev_),
    underRelaxation_(ptf.underRelaxation_),
    aitkenRelaxation_(ptf.aitkenRelaxation_),
    nIter_(ptf.nIter_),
    useAitken_(ptf.useAitken_),
    scalingFactors_(ptf.scalingFactors_),
    overlapPatchNames_(ptf.overlapPatchNames_),
    startTime_(ptf.startTime_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::vector>> Foam::junctionJumpAMIFvPatchVectorField::jump() const
{
    if (this->discontinuousCyclicAMIPatch().owner())
    {
        return jump_;
    }
    else
    {
        const junctionJumpAMIFvPatchVectorField& nbrPatch =
            refCast<const junctionJumpAMIFvPatchVectorField>
            (
                this->neighbourPatchField()
            );

        if (this->discontinuousCyclicAMIPatch().applyLowWeightCorrection())
        {
            return this->discontinuousCyclicAMIPatch().interpolate
            (
                nbrPatch.jump(),
                Field<vector>(this->size(), Zero)
            );
        }
        else
        {

            vectorField normal = this->patch().nf();
            vectorField nbrNormal = this->discontinuousCyclicAMIPatch().neighbPatch().nf();

            // get the jump vector from the neighbour patch
            const vectorField& jmp = nbrPatch.jump();

            // check alignment: sum of (jump · neighbour normal)
            scalar sumDot = sum(jmp & nbrNormal);

            // if anti-parallel, flip the normal
            if (sumDot > 0)
            {
                normal *= -1;
            }

            return this->discontinuousCyclicAMIPatch().interpolate(mag(jmp) * normal);

        }
    }
}


void Foam::junctionJumpAMIFvPatchVectorField::updateCoeffs()
{
    scalar nBranch = overlapPatchNames_.size()+1;

    if(this->discontinuousCyclicAMIPatch().owner())
    {

        
        const customPimpleControl& pimpleLoop = this->db().lookupObject<customPimpleControl>("solutionControl");

        if(pimpleLoop.corr_ ==0) // Check if new timestep
        {
            nIter_ = 0;
        }


        const junctionJumpAMIFvPatchVectorField& nbrPatch =
            refCast<const junctionJumpAMIFvPatchVectorField>
            (
                this->neighbourPatchField()
            );

        // Get patch index from name
        labelList overlapPatchIndeces(nbrPatch.overlappingPatch().size());

        forAll(nbrPatch.overlappingPatch(), i)
        {
            overlapPatchIndeces[i] = this->patch().boundaryMesh().findPatchID(nbrPatch.overlappingPatch()[i]);
            if(overlapPatchIndeces[i] == -1)
            {
                FatalErrorInFunction
                    << "Cannot find patch " << overlapPatchNames_[i] << " in mesh "
                    <<this->patch().boundaryMesh().mesh().name() << exit(FatalError);
            } 
        }

        const surfaceScalarField flux  = this->db().lookupObject<surfaceScalarField>("alphaRhoPhi");

        scalarField restOfFlux(this->size());

        forAll(overlapPatchIndeces, i) // sum of fluxes in all directions except the neighb patch
        {
            restOfFlux += flux.boundaryField()[overlapPatchIndeces[i]];
        }

        const volScalarField& rho=
            this->db().lookupObject<volScalarField>("thermo:rho");

        const scalarField& nbrRho =
            rho.boundaryField()[this->discontinuousCyclicAMIPatch().neighbPatch().index()];
        
        const surfaceScalarField& areas = this->patch().boundaryMesh().mesh().magSf();

        const scalarField& nbrArea =
            areas.boundaryField()[this->discontinuousCyclicAMIPatch().neighbPatch().index()];

        const volScalarField& alpha = 
            this->db().lookupObject<volScalarField>("alpha");

        const scalarField& nbrAlpha =
            alpha.boundaryField()[this->discontinuousCyclicAMIPatch().neighbPatch().index()];   

        scalarField nbrFlux = flux.boundaryField()[this->discontinuousCyclicAMIPatch().neighbPatch().index()];


        // Need to compute the flux like this because if I take the value at the patch
        // it only takes into account the flow through the one patch which is not the total
        // value
        scalarField thisFlux = 
        (
            nBranch*alpha.boundaryField()[this->discontinuousCyclicAMIPatch().index()]
            *rho.boundaryField()[this->discontinuousCyclicAMIPatch().index()]
            *(this->patchInternalField() & this->patch().Sf()) 
        );

        // Compute expected velocity from mass conservation

        scalarField expectedU = -(thisFlux+restOfFlux)/(nbrAlpha*nbrRho*nbrArea+VSMALL); 

        scalarField scalarOwnerU =
        (
            thisFlux[0]>=0 ?
            mag(this->patchInternalField()) :
            -mag(this->patchInternalField())
        );

        scalingFactors_ = nBranch*nbrFlux/(nbrFlux+restOfFlux+VSMALL);
        
        vectorField patchNormal = this->patch().nf();

        vectorField newJump = (-expectedU - scalarOwnerU) * patchNormal;

        // if(useAitken_)
        // {

        //     if (nIter_ == 0)
        //     {
        //         // aitkenRelaxation_ = underRelaxation_;
        //         this->jump_ = this->jump_ + underRelaxation_ * (newJump - this->jump_);
        //         residualPrev_ = newJump[0] - this->jump_[0];
        //     }
        //     else if (nIter_ == 1)
        //     {
        //         this->jump_ = this->jump_ + underRelaxation_ * (newJump - this->jump_);
        //         residual_ = newJump[0] - this->jump_[0];
        //     }
        //     else // nIter_ >= 2 : compute Aitken candidate with robust guards
        //     {
    
        //         vectorField residual(this->size(), residual_);
        //         vectorField residualPrev(this->size(), residualPrev_);
    
        //         scalarField nominator = 
        //         (
        //             residualPrev & (residual - residualPrev)
        //         );
    
        //         scalarField denominator =
        //         (
        //             (residual - residualPrev) & (residual - residualPrev)
        //         );
    
    
        //         aitkenRelaxation_ = mag(-aitkenRelaxation_ * (sum(nominator) / (sum(denominator) + VSMALL)));
    
                
        //         if(aitkenRelaxation_ > 0.05)
        //             aitkenRelaxation_ = 0.05;
            
        //         this->jump_ = this->jump_ + aitkenRelaxation_ * (newJump-this->jump_);
        //         residualPrev_ = residual_;
        //         residual_ = newJump[0] - this->jump_[0];
                
        //     }
    
        //     nIter_++;
    
        // }

        // else
        // {
        this->jump_ = this->jump_ + underRelaxation_ * (newJump - this->jump_);
        // }
    }

    jumpDiscontinuousCyclicAMIFvPatchField<vector>::updateCoeffs();
}



void Foam::junctionJumpAMIFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    jumpDiscontinuousCyclicAMIFvPatchField<vector>::autoMap(m);
    jump_.autoMap(m);
}


void Foam::junctionJumpAMIFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    jumpDiscontinuousCyclicAMIFvPatchField<vector>::rmap(ptf, addr);

    const junctionJumpAMIFvPatchVectorField& tiptf =
        refCast<const junctionJumpAMIFvPatchVectorField>(ptf);
    jump_.rmap(tiptf.jump_, addr);
}


void Foam::junctionJumpAMIFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntry("patchType", this->interfaceFieldType());

    if(this->discontinuousCyclicAMIPatch().owner())
    {
        jump_.writeEntry("jump", os);
    }
    if(!this->discontinuousCyclicAMIPatch().owner())
    {
        os.writeEntry("connectedPatches", overlapPatchNames_);
    }
    scalingFactors_.writeEntry("scalingFactors", os);
    os.writeEntry("underRelaxation", underRelaxation_);
    
    fvPatchField<vector>::writeValueEntry(os);
}



void Foam::junctionJumpAMIFvPatchVectorField::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
 // Collect neighbour face-cell indices
    const labelUList& nbrFaceCells =
        lduAddr.patchAddr
        (
            this->discontinuousCyclicAMIPatch().neighbPatchID()
        );

    // Extract neighbour component field
    solveScalarField pnf(psiInternal, nbrFaceCells);

    // Interpolate across the AMI (handles non-matching meshes)
    pnf = this->discontinuousCyclicAMIPatch().interpolate(pnf);

    // Reference to internal vector field
    const Field<vector>& iField = this->primitiveField();

    // Only apply jump if operating on the original field
    if
    (
        reinterpret_cast<const void*>(&psiInternal)
     == reinterpret_cast<const void*>(&(iField.component(cmpt).ref()))
    )
    {
        Field<vector> jf(this->jump());

        if (!this->discontinuousCyclicAMIPatch().owner())
        {
            jf *= -1.0;
        }

        forAll(pnf, facei)
        {
            pnf[facei] -= jf[facei].component(cmpt);
        }
    }

    // Apply transformation (rotation/reflection) between coupled patches
    this->transformCoupleField(pnf, cmpt);

    // Add into result
    const labelUList& faceCells = lduAddr.patchAddr(patchId);
    this->addToInternalField(result, !add, faceCells, coeffs, pnf);

}


namespace Foam 
{ 
    makePatchTypeField(fvPatchVectorField, junctionJumpAMIFvPatchVectorField); 
}

// ************************************************************************* //
