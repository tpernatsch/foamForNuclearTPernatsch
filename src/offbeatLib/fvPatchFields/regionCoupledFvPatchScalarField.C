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
#include "regionCoupledFvPatchScalarField.H"
#include "Time.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Private members  * * * * * * * * * * * * *//


Foam::tmp<Foam::scalarField>
Foam::regionCoupledFvPatchScalarField::weights() const
{
    if (coupled_)
    {
        const fvPatch& patch = this->patch();

        const scalarField deltas
        (
            patch.nf() & patch.delta()
        );

        const scalarField deltaCoeffs(1/deltas);

        const fvPatch& nbrPatch = regionCoupledPatch_.neighbFvPatch();

        const scalarField nbrDeltas
        (
            regionCoupledPatch_.regionCoupledPatch().interpolate
            (
                nbrPatch.nf() & nbrPatch.delta()
            )
        );

        const scalarField nbrDeltaCoeffs(1/nbrDeltas);

        tmp<scalarField> tw(new scalarField(deltas.size()));
        scalarField& w = tw.ref();

        forAll(deltaCoeffs, facei)
        {
            scalar di = deltaCoeffs[facei];
            scalar dni = nbrDeltaCoeffs[facei];

            w[facei] = di/(di + dni);
        }

        return tw;
    }
    else
    {
        return patch().weights();
    }
}



Foam::tmp<Foam::scalarField>
Foam::regionCoupledFvPatchScalarField::correctionFactors() const
{   
    tmp<Foam::scalarField> tmpCorrectionFactor_(new scalarField(this->size(), 0));
    scalarField& correctionFactor = tmpCorrectionFactor_.ref();

    if (regionCoupledPatch_.owner())
    {

        const AMIInterpolation& ami
            = regionCoupledPatch_.regionCoupledPatch().AMI();

        correctionFactor = ami.srcWeightsSum()*ami.srcMagSf()/patch().magSf();;  

    }
    else
    {
        const AMIInterpolation& ami
            = regionCoupledPatch_.nbrPatch().regionCoupledPatch().AMI();

        correctionFactor = ami.tgtWeightsSum()*ami.tgtMagSf()/patch().magSf();

    }

    return tmpCorrectionFactor_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupledFvPatchScalarField::
regionCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    coupledFvPatchField<scalar>(p, iF),
    regionCoupledPatch_(refCast<const regionCoupledBaseOFFBEATFvPatch>(p)),
    coupled_(false),
    weights_(patch().weights()),
    correctionFactor_(this->size(), 1)
{}


Foam::regionCoupledFvPatchScalarField::
regionCoupledFvPatchScalarField
(
    const regionCoupledFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<scalar>(ptf, p, iF, mapper),
    regionCoupledPatch_(refCast<const regionCoupledBaseOFFBEATFvPatch>(p)),
    coupled_(ptf.coupled_),
    weights_(patch().weights()),
    correctionFactor_(this->size(), 1)
{}


Foam::regionCoupledFvPatchScalarField::
regionCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<scalar>(p, iF, dict),
    regionCoupledPatch_(refCast<const regionCoupledBaseOFFBEATFvPatch>(p)),
    coupled_(dict.lookupOrDefault<bool>("coupled", true)),
    weights_(patch().weights()),
    correctionFactor_(this->size(), 1)
{

    if (!isA<regionCoupledBaseOFFBEAT>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "Foam::regionCoupledFvPatchScalarField::regionCoupledFvPatchScalarField"
            "("
            "    const fvPatch& p,"
            "    const DimensionedField<scalar, volMesh>& iF,"
            "    const dictionary& dict"
            ")"
        )   << "' not type '" << regionCoupledBaseOFFBEAT::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }
}


Foam::regionCoupledFvPatchScalarField::
regionCoupledFvPatchScalarField
(
    const regionCoupledFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    coupledFvPatchField<scalar>(ptf, iF),
    regionCoupledPatch_(ptf.regionCoupledPatch_),
    coupled_(ptf.coupled_),
    weights_(patch().weights()),
    correctionFactor_(this->size(), 1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionCoupledFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    coupledFvPatchField<scalar>::autoMap(m);

#ifdef OPENFOAMFOUNDATION    
    m(weights_, weights_);
    m(correctionFactor_, correctionFactor_);
#elif OPENFOAMESI
    weights_.autoMap(m);
    correctionFactor_.autoMap(m);
#endif    
}


void Foam::regionCoupledFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    coupledFvPatchField<scalar>::rmap(ptf, addr);

    const regionCoupledFvPatchScalarField& dmptf =
        refCast<const regionCoupledFvPatchScalarField>(ptf);

    weights_.rmap(dmptf.weights_, addr);
    correctionFactor_.rmap(dmptf.correctionFactor_, addr);
}


Foam::tmp<Foam::scalarField> Foam::regionCoupledFvPatchScalarField::
snGrad() const
{
    return patch().deltaCoeffs()*(*this - patchInternalField());
}


Foam::tmp<Foam::scalarField> Foam::regionCoupledFvPatchScalarField::
snGrad(const scalarField&) const
{
    return snGrad();
}


void Foam::regionCoupledFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    if (coupled_)
    {
        weights_ = weights();
        correctionFactor_ = correctionFactors();
    }                                                                                                                                                                                                                                                               
    else
    {
        weights_ = patch().weights();

        correctionFactor_ = 1;                                                                                                                                                      
    }
    
    coupledFvPatchField<scalar>::updateCoeffs();
}

void Foam::regionCoupledFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!updated())
    {
        updateCoeffs();
    }

    if (coupled_)
    {
        scalarField::operator=
        (
            weights_*patchInternalField()
            + (1.0 - weights_)*patchNeighbourField()
        );
    }
    else
    {
        // zeroGradient
        scalarField::operator=(patchInternalField());
    }

    fvPatchScalarField::evaluate();
}


Foam::tmp<Foam::Field<Foam::scalar> >
Foam::regionCoupledFvPatchScalarField::
patchNeighbourField() const
{
    const fvPatch& nbrPatch = regionCoupledPatch_.neighbFvPatch();

    const labelUList& nbrFaceCells = nbrPatch.faceCells();
    
    const scalarField npf(internalField(), nbrFaceCells);   

    return regionCoupledPatch_.regionCoupledPatch().interpolate(npf);
}


void Foam::regionCoupledFvPatchScalarField::
patchNeighbourField(UList<scalar>& pnf) const
{
    const Field<scalar>& iField = this->primitiveField();
    const labelUList& nbrFaceCells =
        regionCoupledPatch_.neighbFvPatch().faceCells();

    forAll(pnf, facei)
    {
        pnf[facei] = iField[nbrFaceCells[facei]];
    }
}


Foam::tmp<Foam::scalarField>
Foam::regionCoupledFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    if (coupled_)
    {
        return weights_;
    }
    else
    {
        return tmp<scalarField>
        (
            new scalarField(this->size(), 1)
        );
    }
}


Foam::tmp<Foam::scalarField>
Foam::regionCoupledFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    if (coupled_)
    {
        return (1.0 - weights_);
    }
    else
    {
        return tmp<scalarField>
        (
            new scalarField(this->size(), 0)
        );
    }
}


Foam::tmp<Foam::scalarField>
Foam::regionCoupledFvPatchScalarField::gradientInternalCoeffs
(
    const scalarField&
) const
{
    return gradientInternalCoeffs();
}


Foam::tmp<Foam::scalarField>
Foam::regionCoupledFvPatchScalarField::gradientInternalCoeffs() const
{
    if (coupled_)
    {
        return correctionFactor_*(weights_ - 1)
            *regionCoupledPatch_.patch().deltaCoeffs();
    }
    else
    {
        return tmp<scalarField>
        (
            new scalarField(this->size(), 0)
        );
    }
}


Foam::tmp<Foam::scalarField>
Foam::regionCoupledFvPatchScalarField::gradientBoundaryCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    return gradientBoundaryCoeffs();
}


Foam::tmp<Foam::scalarField>
Foam::regionCoupledFvPatchScalarField::gradientBoundaryCoeffs() const
{
    if (coupled_)
    {
        return correctionFactor_*(1 - weights_)
            *regionCoupledPatch_.patch().deltaCoeffs();
    }
    else
    {
        return tmp<scalarField>
        (
            new scalarField(this->size(), 0)
        );
    }
}

#ifdef OPENFOAMFOUNDATION
void Foam::regionCoupledFvPatchScalarField::updateInterfaceMatrix
(
    Field<scalar>& result,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes comms
) const
{
    updateInterfaceMatrix(result, psiInternal, coeffs, comms);
}


void Foam::regionCoupledFvPatchScalarField::updateInterfaceMatrix
(
    Field<scalar>& result,
    const Field<scalar>& psiInternal,
    const scalarField& coeffs,
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

#elif OPENFOAMESI
void Foam::regionCoupledFvPatchScalarField::updateInterfaceMatrix
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


void Foam::regionCoupledFvPatchScalarField::updateInterfaceMatrix
(
Field<scalar>& result,
const bool add,
const lduAddressing& lduAddr,
const label patchId,
const Field<scalar>& psiInternal,
const scalarField& coeffs,
const Pstream::commsTypes
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
    this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}
#endif


void Foam::regionCoupledFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
#ifdef OPENFOAMFOUNDATION    
    writeEntry(os, "patchType", regionCoupledPatch_.type());
    writeEntry(os, "coupled", coupled_);
    writeEntry(os, "value", *this);
#elif OPENFOAMESI
    os.writeEntry("patchType", regionCoupledPatch_.type());
    os.writeEntry("coupled", coupled_);
    this->writeEntry("value", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        regionCoupledFvPatchScalarField
    );
};


// ************************************************************************* //
