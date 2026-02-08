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

#include "resistiveLayerFvPatchScalarField.H"
#include "scalarFieldFieldINew.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{
 
//* * * * * * * * * * * * * Protected member functions  * * * * * * * * * * *//

void resistiveLayerFvPatchScalarField::check() const
{
    dimensionSet dimConductivity = dimPower/dimLength/dimTemperature;
    
    if (internalField().dimensions() != dimConductivity)
    {
        FatalErrorInFunction()
            << "The resistiveLayer boundary is intended for use "
               "on a thermal conductivity field." << nl
            << "Field " << internalField().name() << " dimensions: "
            << internalField().dimensions() << nl
            << "Expected dimensions: " << dimConductivity << nl
            << abort(FatalError);
    }
}

void resistiveLayerFvPatchScalarField::applyCorrection()
{
    const scalarField rDelta = patch().deltaCoeffs();
    scalarField alpha_bulk = (*this)*rDelta;
    scalarField alpha_eff = alpha_bulk / (1 + alpha_bulk/max(alpha_, SMALL));

    fvPatchScalarField::operator=(alpha_eff / rDelta);
}


//  * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //
 
resistiveLayerFvPatchScalarField::resistiveLayerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    alpha_(p.size(), GREAT)
{
    check();
}

resistiveLayerFvPatchScalarField::resistiveLayerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    alpha_(p.size(), GREAT)
{
    check();
    
    if(dict.found("alpha"))
    {
        alpha_ = scalarField("alpha", dict, p.size());
    }
}

 
resistiveLayerFvPatchScalarField::resistiveLayerFvPatchScalarField
(
    const resistiveLayerFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    alpha_(mapper(ptf.alpha_))
{
    check();
}
 
 
resistiveLayerFvPatchScalarField::
resistiveLayerFvPatchScalarField
(
    const resistiveLayerFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    alpha_(ptf.alpha_)
{
    check();
}
 

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void resistiveLayerFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
#ifdef OPENFOAMFOUNDATION
    m(alpha_, alpha_);
#elif OPENFOAMESI
    alpha_.autoMap(m);
#endif  
}


void resistiveLayerFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const resistiveLayerFvPatchScalarField& mptf =
     refCast<const resistiveLayerFvPatchScalarField>(ptf);

    alpha_.rmap(mptf.alpha_, addr);
}


void resistiveLayerFvPatchScalarField::updateCoeffs()
{
    this->operator=(patchInternalField());
}


void resistiveLayerFvPatchScalarField::operator=(const UList<scalar>& pf)
{
    fvPatchScalarField::operator=(pf);
    this->applyCorrection();
}


void resistiveLayerFvPatchScalarField::operator=(const fvPatchScalarField& pf)
{
    fvPatchScalarField::operator=(pf);
    this->applyCorrection();
}


void resistiveLayerFvPatchScalarField::operator=(const scalar& value)
{
    fvPatchScalarField::operator=(value);
    this->applyCorrection();
}


void resistiveLayerFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);

    #ifdef OPENFOAMFOUNDATION
    writeEntry(os, "alpha", alpha_);
    
    #elif OPENFOAMESI
    alpha_.writeEntry("alpha", os);
    
    #endif    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    resistiveLayerFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam