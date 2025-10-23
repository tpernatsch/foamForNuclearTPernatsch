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

#include "branchingJunctionJumpAMIFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::branchingJunctionJumpAMIFvPatchVectorField::branchingJunctionJumpAMIFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    jumpDiscontinuousCyclicAMIFvPatchField<vector>(p, iF),
    jump_(this->size(), Zero),
    underRelaxation_(1),
    scalingFactors_(this->size(), Zero),
    overlapPatchNames_(this->size(), "")
{}


Foam::branchingJunctionJumpAMIFvPatchVectorField::branchingJunctionJumpAMIFvPatchVectorField
(
    const branchingJunctionJumpAMIFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    jumpDiscontinuousCyclicAMIFvPatchField<vector>(ptf, p, iF, mapper),
    jump_(ptf.jump_, mapper),
    underRelaxation_(ptf.underRelaxation_),
    scalingFactors_(ptf.scalingFactors_),
    overlapPatchNames_(ptf.overlapPatchNames_)
{}


Foam::branchingJunctionJumpAMIFvPatchVectorField::branchingJunctionJumpAMIFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    jumpDiscontinuousCyclicAMIFvPatchField<vector>(p, iF),
    jump_(p.size(), Zero),
    underRelaxation_(dict.getOrDefault<scalar>("underRelaxation", 1)),
    scalingFactors_(this->size(), Zero),
    overlapPatchNames_(this->size(), "")
{
    if(!this->discontinuousCyclicAMIPatch().owner())
    {
        overlapPatchNames_ = dict.get<wordList>("connectedPatches");
    }
}



Foam::branchingJunctionJumpAMIFvPatchVectorField::branchingJunctionJumpAMIFvPatchVectorField
(
    const branchingJunctionJumpAMIFvPatchVectorField& ptf
)
:
    jumpDiscontinuousCyclicAMIFvPatchField<vector>(ptf),
    jump_(ptf.jump_),
    underRelaxation_(ptf.underRelaxation_),
    scalingFactors_(ptf.scalingFactors_),
    overlapPatchNames_(ptf.overlapPatchNames_)
{}


Foam::branchingJunctionJumpAMIFvPatchVectorField::branchingJunctionJumpAMIFvPatchVectorField
(
    const branchingJunctionJumpAMIFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    jumpDiscontinuousCyclicAMIFvPatchField<vector>(ptf, iF),
    jump_(ptf.jump_),
    underRelaxation_(ptf.underRelaxation_),
    scalingFactors_(ptf.scalingFactors_),
    overlapPatchNames_(ptf.overlapPatchNames_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::vector>> Foam::branchingJunctionJumpAMIFvPatchVectorField::jump() const
{
    if (this->discontinuousCyclicAMIPatch().owner())
    {
        return jump_;
    }
    else
    {
        const branchingJunctionJumpAMIFvPatchVectorField& nbrPatch =
            refCast<const branchingJunctionJumpAMIFvPatchVectorField>
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


void Foam::branchingJunctionJumpAMIFvPatchVectorField::updateCoeffs()
{
    if(this->discontinuousCyclicAMIPatch().owner())
    {
        const branchingJunctionJumpAMIFvPatchVectorField& nbrPatch =
            refCast<const branchingJunctionJumpAMIFvPatchVectorField>
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
            restOfFlux += mag(flux.boundaryField()[overlapPatchIndeces[i]]);
        }

        scalarField thisFlux = mag(flux.boundaryField()[this->discontinuousCyclicAMIPatch().index()]); 
        scalarField nbrFlux = mag(flux.boundaryField()[this->discontinuousCyclicAMIPatch().neighbPatch().index()]);

        scalarField nominator = nbrFlux;
        scalarField denominator = restOfFlux + nbrFlux + VSMALL;

        // Same for density and areas

        const volScalarField& rho=
            this->db().lookupObject<volScalarField>("thermo:rho");

        const scalarField& nbrRho =
            rho.boundaryField()[this->discontinuousCyclicAMIPatch().neighbPatch().index()];
        const scalarField& thisRho =
            rho.boundaryField()[this->discontinuousCyclicAMIPatch().index()];

        
        const surfaceScalarField& areas = this->patch().boundaryMesh().mesh().magSf();

        const scalarField& nbrArea =
            areas.boundaryField()[this->discontinuousCyclicAMIPatch().neighbPatch().index()];
        const scalarField& thisArea =
            areas.boundaryField()[this->discontinuousCyclicAMIPatch().index()];

        
        const volScalarField& alpha = 
            this->db().lookupObject<volScalarField>("alpha");

        const scalarField& nbrAlpha =
            alpha.boundaryField()[this->discontinuousCyclicAMIPatch().neighbPatch().index()];   
        const scalarField& thisAlpha =
            alpha.boundaryField()[this->discontinuousCyclicAMIPatch().index()];


        if(min(denominator) == 0)
        {
  
            this->jump_ = 0*this->patchInternalField();

        }
        else
        {
            scalingFactors_ = nominator/denominator;    

            vectorField velocityVersor = this->patchInternalField()/(mag(this->patchInternalField())+VSMALL);

            vectorField LHS = this->patchInternalField()*(thisRho*thisAlpha*thisArea- nbrRho*nbrArea*nbrAlpha)/(nbrRho*nbrArea*nbrAlpha+VSMALL);
            vectorField RHS = restOfFlux*velocityVersor/(nbrRho*nbrArea*nbrAlpha+VSMALL);
        
            this->jump_ = 
            (
                underRelaxation_*(LHS-RHS) + (1-underRelaxation_)*this->jump_   

            );
        }

    }

    jumpDiscontinuousCyclicAMIFvPatchField<vector>::updateCoeffs();
}



void Foam::branchingJunctionJumpAMIFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    jumpDiscontinuousCyclicAMIFvPatchField<vector>::autoMap(m);
    jump_.autoMap(m);
}


void Foam::branchingJunctionJumpAMIFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    jumpDiscontinuousCyclicAMIFvPatchField<vector>::rmap(ptf, addr);

    const branchingJunctionJumpAMIFvPatchVectorField& tiptf =
        refCast<const branchingJunctionJumpAMIFvPatchVectorField>(ptf);
    jump_.rmap(tiptf.jump_, addr);
}


void Foam::branchingJunctionJumpAMIFvPatchVectorField::write(Ostream& os) const
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
    
    fvPatchField<vector>::writeValueEntry(os);
}



void Foam::branchingJunctionJumpAMIFvPatchVectorField::updateInterfaceMatrix
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
    makePatchTypeField(fvPatchVectorField, branchingJunctionJumpAMIFvPatchVectorField); 
}

// ************************************************************************* //
