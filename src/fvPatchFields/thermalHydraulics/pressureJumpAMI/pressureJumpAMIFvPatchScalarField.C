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

#include "pressureJumpAMIFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureJumpAMIFvPatchScalarField::pressureJumpAMIFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    jumpDiscontinuousCyclicAMIFvPatchField<scalar>(p, iF),
    jump_(this->size(), Zero),
    bernoulliCorrection_(false),
    source_(this->size(), Zero),
    sourceTable_(nullptr),
    pressureLossCoeff_(0),
    pumpCoeffs_(1, Zero),
    pumpVelocityTable_(nullptr),
    nominalVelocity_(1)
{}


Foam::pressureJumpAMIFvPatchScalarField::pressureJumpAMIFvPatchScalarField
(
    const pressureJumpAMIFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    jumpDiscontinuousCyclicAMIFvPatchField<scalar>(ptf, p, iF, mapper),
    jump_(ptf.jump_, mapper),
    bernoulliCorrection_(ptf.bernoulliCorrection_),
    source_(ptf.source_, mapper),
    sourceTable_(ptf.sourceTable_.clone()),
    pressureLossCoeff_(ptf.pressureLossCoeff_),
    pumpCoeffs_(ptf.pumpCoeffs_),
    pumpVelocityTable_(ptf.pumpVelocityTable_.clone()),
    nominalVelocity_(ptf.nominalVelocity_)
{}


Foam::pressureJumpAMIFvPatchScalarField::pressureJumpAMIFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    jumpDiscontinuousCyclicAMIFvPatchField<scalar>(p, iF),
    jump_(p.size(), Zero),
    bernoulliCorrection_(dict.getOrDefault<bool>("bernoulli", true)),
    source_(p.size(), Zero),
    sourceTable_(nullptr),
    pressureLossCoeff_(dict.getOrDefault<scalar>("lossCoeff", 0)),
    pumpCoeffs_(dict.getOrDefault<scalarList>("pumpCoeffs", scalarList(1, Zero))),
    pumpVelocityTable_(nullptr),
    nominalVelocity_(1)
{
    source_.assign("pSource", dict, p.size(), IOobjectOption::LAZY_READ);

    if(dict.found("pumpVelocityTable"))
    {
        pumpVelocityTable_ = Function1<scalar>::New("pumpVelocityTable", dict, &this->db());
        nominalVelocity_ = dict.get<scalar>("nominalVelocity");
    }

    if(dict.found("sourceTable"))
    {
        sourceTable_ = Function1<scalar>::New("sourceTable", dict, &this->db());
    }


}


Foam::pressureJumpAMIFvPatchScalarField::pressureJumpAMIFvPatchScalarField
(
    const pressureJumpAMIFvPatchScalarField& ptf
)
:
    jumpDiscontinuousCyclicAMIFvPatchField<scalar>(ptf),
    jump_(ptf.jump_),
    bernoulliCorrection_(ptf.bernoulliCorrection_),
    source_(ptf.source_),
    sourceTable_(ptf.sourceTable_.clone()),
    pressureLossCoeff_(ptf.pressureLossCoeff_),
    pumpCoeffs_(ptf.pumpCoeffs_),
    pumpVelocityTable_(ptf.pumpVelocityTable_.clone()),
    nominalVelocity_(ptf.nominalVelocity_)
{}


Foam::pressureJumpAMIFvPatchScalarField::pressureJumpAMIFvPatchScalarField
(
    const pressureJumpAMIFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    jumpDiscontinuousCyclicAMIFvPatchField<scalar>(ptf, iF),
    jump_(ptf.jump_),
    bernoulliCorrection_(ptf.bernoulliCorrection_),
    source_(ptf.source_),
    sourceTable_(ptf.sourceTable_.clone()),
    pressureLossCoeff_(ptf.pressureLossCoeff_),
    pumpCoeffs_(ptf.pumpCoeffs_),
    pumpVelocityTable_(ptf.pumpVelocityTable_.clone()),
    nominalVelocity_(ptf.nominalVelocity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::scalar>> Foam::pressureJumpAMIFvPatchScalarField::jump() const
{
    if (this->discontinuousCyclicAMIPatch().owner())
    {
        return jump_;
    }
    else
    {

        const pressureJumpAMIFvPatchScalarField& nbrPatch =
            refCast<const pressureJumpAMIFvPatchScalarField>
            (
                this->neighbourPatchField()
            );

        if (this->discontinuousCyclicAMIPatch().applyLowWeightCorrection())
        {
            return this->discontinuousCyclicAMIPatch().interpolate
            (
                nbrPatch.jump(),
                Field<scalar>(this->size(), Zero)
            );
        }
        else
        {
            return this->discontinuousCyclicAMIPatch().interpolate(nbrPatch.jump());
        }
    }
}

void Foam::pressureJumpAMIFvPatchScalarField::updateCoeffs()
{

    if (this->updated())
    {
        return;
    }

    if(this->discontinuousCyclicAMIPatch().owner())
    {

        // Update jump to account for differences in velocity

        // In order to account for the correct flow direction,
        // first I check the flux value. If positive the flow is 
        // from the owner into the slave, hence vOwner is employed
        // for the calculations. Otherwise, I have to do the opposite

        const scalarField& fluxOwner = 
            this->patch().lookupPatchField<surfaceScalarField, scalar>("alphaRhoPhi");

        bool ownerToSlave = (fluxOwner[0]>=0);

        const vectorField& vOwner  = 
            this->patch().lookupPatchField<volVectorField, vector>("U");

        const scalarField& rhoOwner=
            this->patch().lookupPatchField<volScalarField, scalar>("thermo:rho");

        const scalarField& rhoSlave=
            this->discontinuousCyclicAMIPatch().neighbPatch().lookupPatchField<volScalarField, scalar>("thermo:rho");

        scalarField rhoSlaveOnMaster(vOwner.size());

        rhoSlaveOnMaster = this->discontinuousCyclicAMIPatch().interpolate(rhoSlave);

        const vectorField& vSlave = 
            this->discontinuousCyclicAMIPatch().neighbPatch().lookupPatchField<volVectorField, vector>("U");

        scalarField magUOwner = mag(vOwner);
        scalarField magUSlave = mag(vSlave);

        magUSlave = this->discontinuousCyclicAMIPatch().interpolate(magUSlave); 

        // Get massflow rate on master patch

        scalarField pumpSource(vOwner.size(), Zero);

        if(pumpCoeffs_.size() > 1 or pumpCoeffs_[0] != 0)
        {

            scalarField massFlowOwner
                = this->patch().lookupPatchField<surfaceScalarField, scalar>("alphaPhi");

            scalarList correctedCoeffs (pumpCoeffs_);

            if(pumpVelocityTable_ != nullptr)
            {
                scalar newSpeed = pumpVelocityTable_->value(this->db().time().value());
                scalar speedRatio = newSpeed/nominalVelocity_;

                if(speedRatio < SMALL)
                {
                    correctedCoeffs = scalarList(pumpCoeffs_.size(), 0);
                }

                else
                {
                    forAll(correctedCoeffs, coeffI)
                    {
                        correctedCoeffs[coeffI] = pumpCoeffs_[coeffI]*pow(speedRatio, 2-coeffI);
                    }
                }
            }

            uniformDimensionedVectorField g = this->db().lookupObject<uniformDimensionedVectorField>("g");
            

            forAll(correctedCoeffs, coeffI)
            {
                pumpSource += correctedCoeffs[coeffI] * pow(massFlowOwner, coeffI) * rhoOwner * mag(g.value());
            }
        }

        scalarField correctedSource(source_);

        if(sourceTable_ != nullptr)
        {
            
            correctedSource = source_ * sourceTable_->value(this->db().time().value());
            
        }

        this->jump_ = 
        (
            ownerToSlave? 
                correctedSource - pressureLossCoeff_*0.5*rhoOwner*pow(magUOwner, 2) + pumpSource :
                correctedSource + pressureLossCoeff_*0.5*rhoSlaveOnMaster*pow(magUSlave,2) + pumpSource
        );

        if(bernoulliCorrection_)
        {  
            this->jump_ +=0.5*(rhoOwner*pow(magUOwner,2) - rhoSlaveOnMaster*pow(magUSlave,2)); // bernoulli
        }

    }

    jumpDiscontinuousCyclicAMIFvPatchField<scalar>::updateCoeffs();

}



void Foam::pressureJumpAMIFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    jumpDiscontinuousCyclicAMIFvPatchField<scalar>::autoMap(m);
    jump_.autoMap(m);
}


void Foam::pressureJumpAMIFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    jumpDiscontinuousCyclicAMIFvPatchField<scalar>::rmap(ptf, addr);

    const pressureJumpAMIFvPatchScalarField& tiptf =
        refCast<const pressureJumpAMIFvPatchScalarField>(ptf);
    jump_.rmap(tiptf.jump_, addr);
}


void Foam::pressureJumpAMIFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("patchType", this->interfaceFieldType());

    if (this->discontinuousCyclicAMIPatch().owner())
    {
        jump_.writeEntry("jump", os);
        source_.writeEntry("pSource", os);
    }

    fvPatchField<scalar>::writeValueEntry(os);
}


namespace Foam 
{ 
    makePatchTypeField(fvPatchScalarField, pressureJumpAMIFvPatchScalarField); 
}

// ************************************************************************* //
