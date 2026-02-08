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

#include "gaussianLaserFvPatchScalarField.H"
#include "scalarFieldFieldINew.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{
 
//* * * * * * * * * * * * * Protected member functions  * * * * * * * * * * *//
 
//  * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //
 
gaussianLaserFvPatchScalarField::gaussianLaserFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    P_(0.0),
    w_(0.0),
    R_(0.0),
    powerList_()
{}

gaussianLaserFvPatchScalarField::gaussianLaserFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict, false),
    P_(0.0),
    w_(readScalar(dict.lookup("beamWaist"))),
    R_(readScalar(dict.lookup("reflectionCoeff"))),
    powerList_()
{
    // Read LASER power (either list or fixed)
    if(dict.found("powerList"))
    {
#ifdef OPENFOAMFOUNDATION            
        powerList_.set
        ( 
            new Function1s::Table<scalar>
            (
                "powerList", dict.subDict("powerList")
            )
        );
#elif OPENFOAMESI            
        powerList_.reset
        ( 
            new scalarTable
            (
                "powerList", dict.subDict("powerList")
            )
        );
#endif                        

        P_ = 
        powerList_->value
        (
            patch().boundaryMesh().mesh().time().timeToUserTime
            (
                this->db().time().value()
            )
        );    
    }
    else if (dict.found("power"))
    {
        P_ = readScalar(dict.lookup("power"));
    }
    else
    {
        FatalErrorInFunction() 
        << "gaussianLaser BC for " << patch().name() << " requires:" << nl
        << "- \"powerList\" or" << nl
        << "- \"power\"" << abort(FatalError) << endl;
    }
}
 
 
gaussianLaserFvPatchScalarField::gaussianLaserFvPatchScalarField
(
    const gaussianLaserFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
#ifdef OPENFOAMFOUNDATION
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper, mappingRequired),
#elif OPENFOAMESI
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
#endif    
    P_(ptf.P_),
    w_(ptf.w_),
    R_(ptf.R_)
{}
 
 
gaussianLaserFvPatchScalarField::
gaussianLaserFvPatchScalarField
(
    const gaussianLaserFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    P_(ptf.P_),
    w_(ptf.w_),
    R_(ptf.R_),
    powerList_()
{}
 
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
void Foam::gaussianLaserFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if(powerList_.valid())
    {
        P_ = 
        powerList_->value
        (
            patch().boundaryMesh().mesh().time().timeToUserTime
            (
                this->db().time().value()
            )
        );
    }

    // Initialize patch field
    scalarField Ip(this->size(), 0.0);

    // Compute radius of patch faces (assumed r-theta geometry)
    scalarField r(this->size(), 0.0);
    forAll(patch().Cf(), i)
    {
        r[i] = pow
        (
            pow(patch().Cf()[i].x(), 2.0) + pow(patch().Cf()[i].y(), 2.0)
            , 
            0.5
        );
    }

    // Compute intensity at this patch
    Ip = 2*P_/(3.1415*pow(w_, 2.0)) * exp(-2*pow(r,2.0)/pow(w_, 2.0)) * (1-R_);

    // Assign to patch field
    fvPatchField<scalar>::operator==(Ip);

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::gaussianLaserFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

#ifdef OPENFOAMFOUNDATION
    writeEntry(os, "value", *this);
    writeEntry(os, "beamWaist", w_);
    writeEntry(os, "power", P_);
    writeEntry(os, "reflectionCoeff", R_);
#elif OPENFOAMESI
    this->writeEntry("value", os);    
    os.writeEntry("beamWaist", w_);
    os.writeEntry("power", P_);
    os.writeEntry("reflectionCoeff", R_);
#endif    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    gaussianLaserFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam