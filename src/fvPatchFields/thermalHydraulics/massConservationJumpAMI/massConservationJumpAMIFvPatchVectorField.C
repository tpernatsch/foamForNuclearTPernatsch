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

#include "massConservationJumpAMIFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massConservationJumpAMIFvPatchVectorField::massConservationJumpAMIFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    jumpDiscontinuousCyclicAMIFvPatchField<vector>(p, iF),
    jump_(this->size(), Zero),
    underRelaxation_(1),
    scalingFactors_(this->size(), Zero),
    startTime_(0)
{}


Foam::massConservationJumpAMIFvPatchVectorField::massConservationJumpAMIFvPatchVectorField
(
    const massConservationJumpAMIFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    jumpDiscontinuousCyclicAMIFvPatchField<vector>(ptf, p, iF, mapper),
    jump_(ptf.jump_, mapper),
    underRelaxation_(ptf.underRelaxation_),
    scalingFactors_(ptf.scalingFactors_),
    startTime_(ptf.startTime_)
{}


Foam::massConservationJumpAMIFvPatchVectorField::massConservationJumpAMIFvPatchVectorField
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
    startTime_(dict.getOrDefault<scalar>("startTime", 0))

{
    if (this->discontinuousCyclicAMIPatch().owner())
    {
        jump_.assign("jump", dict, p.size(), IOobjectOption::LAZY_READ);
    }
}


Foam::massConservationJumpAMIFvPatchVectorField::massConservationJumpAMIFvPatchVectorField
(
    const massConservationJumpAMIFvPatchVectorField& ptf
)
:
    jumpDiscontinuousCyclicAMIFvPatchField<vector>(ptf),
    jump_(ptf.jump_),
    underRelaxation_(ptf.underRelaxation_),
    scalingFactors_(ptf.scalingFactors_),
    startTime_(ptf.startTime_)
{}


Foam::massConservationJumpAMIFvPatchVectorField::massConservationJumpAMIFvPatchVectorField
(
    const massConservationJumpAMIFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    jumpDiscontinuousCyclicAMIFvPatchField<vector>(ptf, iF),
    jump_(ptf.jump_),
    underRelaxation_(ptf.underRelaxation_),
    scalingFactors_(ptf.scalingFactors_),
    startTime_(ptf.startTime_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::vector>> Foam::massConservationJumpAMIFvPatchVectorField::jump() const
{
    if (this->discontinuousCyclicAMIPatch().owner())
    {
        return jump_;
    }
    else
    {
        const massConservationJumpAMIFvPatchVectorField& nbrPatch =
            refCast<const massConservationJumpAMIFvPatchVectorField>
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


// void Foam::massConservationJumpAMIFvPatchVectorField::evaluate(const Pstream::commsTypes)
// {


//     if (!this->discontinuousCyclicAMIPatch().owner())
//     {
    
//         fvPatchField<vector>::operator == 
//         (
//             this->discontinuousCyclicAMIPatch().neighbPatch().lookupPatchField<volVectorField, vector>("U")
//         );
//     }

// }

void Foam::massConservationJumpAMIFvPatchVectorField::updateCoeffs()
{

    if(this->db().time().value()>=startTime_)
    {
        if(this->discontinuousCyclicAMIPatch().owner())
        {
            // Update jump to account for differences in velocity
    
            const vectorField& vOwner = 
                this->patch().lookupPatchField<volVectorField, vector>("U");
    
            const scalarField& rhoOwner=
                this->patch().lookupPatchField<volScalarField, scalar>("thermo:rho");
    
            const scalarField& rhoSlave =
                this->discontinuousCyclicAMIPatch().neighbPatch().lookupPatchField<volScalarField, scalar>("thermo:rho");
            scalarField rhoSlaveOnMaster(vOwner.size());
            rhoSlaveOnMaster = this->discontinuousCyclicAMIPatch().interpolate(rhoSlave);
    
            const scalarField& ownerAreas = this->patch().magSf();
    
            const scalarField slaveAreas = this->discontinuousCyclicAMIPatch().neighbPatch().magSf();
            scalarField slaveAreasSumField(slaveAreas.size(), gSum(slaveAreas));
            scalarField ownerAreasSumField(ownerAreas.size(), gSum(ownerAreas));
    
            scalarField slaveAreasOnOwner(vOwner.size());
            slaveAreasOnOwner = this->discontinuousCyclicAMIPatch().interpolate(slaveAreasSumField);
    
            scalingFactors_ = rhoOwner*ownerAreasSumField/rhoSlaveOnMaster/slaveAreasOnOwner;
    
    
    
    
            // Adjust for porosity
    
            volScalarField alpha =
                this->db().lookupObject<volScalarField>("alpha");
            
            scalarField alphaOwner =
                alpha.boundaryField()[ this->patch().index() ].patchInternalField();
    
            scalarField alphaSlave = 
                alpha.boundaryField()[ this->discontinuousCyclicAMIPatch().neighbPatch().index() ].patchInternalField();
    
            scalarField alphaSlaveOnMaster(vOwner.size());
    
            alphaSlaveOnMaster = this->discontinuousCyclicAMIPatch().interpolate(alphaSlave);
    
    
            scalingFactors_ = scalingFactors_ * (alphaOwner/alphaSlaveOnMaster);
    
    
            this->jump_ = 
                underRelaxation_*this->patchInternalField()*(scalingFactors_ -1)
               +(1-underRelaxation_)*this->jump_;          
            // this->jump_ = this->patchInternalField()*0;
    
    
        }
    
        jumpDiscontinuousCyclicAMIFvPatchField<vector>::updateCoeffs();
    }

}



void Foam::massConservationJumpAMIFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    jumpDiscontinuousCyclicAMIFvPatchField<vector>::autoMap(m);
    jump_.autoMap(m);
}


void Foam::massConservationJumpAMIFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    jumpDiscontinuousCyclicAMIFvPatchField<vector>::rmap(ptf, addr);

    const massConservationJumpAMIFvPatchVectorField& tiptf =
        refCast<const massConservationJumpAMIFvPatchVectorField>(ptf);
    jump_.rmap(tiptf.jump_, addr);
}


void Foam::massConservationJumpAMIFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntry("patchType", this->interfaceFieldType());

   
    jump_.writeEntry("jump", os);
    scalingFactors_.writeEntry("scalingFactors", os);

    os.writeEntry("underRelaxation", underRelaxation_);
    
    fvPatchField<vector>::writeValueEntry(os);
}



void Foam::massConservationJumpAMIFvPatchVectorField::updateInterfaceMatrix
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




// Foam::tmp<Foam::Field<Foam::vector>>
// Foam::massConservationJumpAMIFvPatchVectorField::valueInternalCoeffs
// (
//     const tmp<scalarField>& w
// ) const
// {
//     Foam::vector val = this->discontinuousCyclicAMIPatch().owner()
//         ? Foam::pTraits<Foam::vector>::one
//         : Foam::pTraits<Foam::vector>::zero;

//     return Foam::tmp<Foam::Field<Foam::vector>>
//     (
//         new Foam::Field<Foam::vector>(this->patch().size(), val)
//     );
// }


// Foam::tmp<Foam::Field<Foam::vector>>
// Foam::massConservationJumpAMIFvPatchVectorField::valueBoundaryCoeffs
// (
//     const tmp<scalarField>& w
// ) const
// {

//     Foam::vector val = this->discontinuousCyclicAMIPatch().owner()
//         ? Foam::pTraits<Foam::vector>::zero
//         : Foam::pTraits<Foam::vector>::one;

//     return Foam::tmp<Foam::Field<Foam::vector>>
//     (
//         new Foam::Field<Foam::vector>(this->patch().size(), val)
//     );
// }


// Foam::tmp<Foam::Field<Foam::vector>>
// Foam::massConservationJumpAMIFvPatchVectorField::gradientInternalCoeffs
// (
//     const scalarField& deltaCoeffs
// ) const
// {
//     return -vector(pTraits<vector>::one)*deltaCoeffs;
// }


// Foam::tmp<Foam::Field<Foam::vector>>
// Foam::massConservationJumpAMIFvPatchVectorField::gradientInternalCoeffs() const
// {
//     return -vector(pTraits<vector>::one)*this->patch().deltaCoeffs();
// }


// Foam::tmp<Foam::Field<Foam::vector>>
// Foam::massConservationJumpAMIFvPatchVectorField::gradientBoundaryCoeffs
// (
//     const scalarField& deltaCoeffs
// ) const
// {
//     return -this->gradientInternalCoeffs(deltaCoeffs);
// }


// Foam::tmp<Foam::Field<Foam::vector>>
// Foam::massConservationJumpAMIFvPatchVectorField::gradientBoundaryCoeffs() const
// {
//     return -this->gradientInternalCoeffs();
// }





// void Foam::massConservationJumpAMIFvPatchVectorField::evaluate(const Pstream::commsTypes) 
// {
//     if (!this->updated())
//     {
//         this->updateCoeffs();
//     }

//     if(!this->cyclicAMIFvPatch().owner())
//     {
//         // Get other velocity value

//         const vectorField& vOwner = 
//             this->patch().lookupPatchField<volVectorField, vector>("U");
//     }


// }






namespace Foam 
{ 
    makePatchTypeField(fvPatchVectorField, massConservationJumpAMIFvPatchVectorField); 
}

// ************************************************************************* //
