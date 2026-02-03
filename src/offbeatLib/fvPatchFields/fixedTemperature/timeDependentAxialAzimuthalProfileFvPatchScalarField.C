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

#include "timeDependentAxialAzimuthalProfileFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "scalarFieldFieldINew.H"

namespace Foam
{

//  * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void timeDependentAxialAzimuthalProfileFvPatchScalarField::calculateTheta()
{
    const vectorField& Cf = referencePatch_.Cf();
    rAxis_ /= mag(rAxis_);
    
    vectorField R = Cf - z_*axis_;
    
    theta_ = atan2(axis_ & (rAxis_ ^ R), rAxis_ & R)
             + constant::mathematical::pi;
}

//  * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
timeDependentAxialAzimuthalProfileFvPatchScalarField::
timeDependentAxialAzimuthalProfileFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    timeDependentAxialProfileFvPatchScalarField(p, iF),
    thetaValues_(),
	thetaData_(),
    thetaMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    thetaTable_(timeValues_, thetaData_, thetaMethod_),
    rAxis_(1,0,0),
    theta_()
{
    if (db().foundObject<globalOptions>("globalOptions"))
    {
        rAxis_ = db().lookupObject<globalOptions>("globalOptions").radialDirection();
    }

    calculateTheta();
}

timeDependentAxialAzimuthalProfileFvPatchScalarField::
timeDependentAxialAzimuthalProfileFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    timeDependentAxialProfileFvPatchScalarField(p, iF, dict, valueRequired),
    thetaValues_(profileDict_.lookup("thetaLocations")),
    thetaData_
    (
        PtrList<scalarField>
        (
            profileDict_.lookup("thetaData"), 
            scalarFieldFieldINew()
        )
    ),
    thetaMethod_(interpolateTableBase::interpolationMethodNames_[
        profileDict_.lookupOrDefault<word>("thetaInterpolationMethod", "linear")
    ]),
    thetaTable_(timeValues_, thetaData_, thetaMethod_),
    rAxis_(1,0,0),
    theta_()
{
    if ((min(thetaValues_) < 0) or (max(thetaValues_) > 360))
    {
        FatalErrorInFunction()
            << "Supplied angular locations in thetaLocations must lie in "
            << "the range [0, 360]." << endl
            << abort(FatalError);
    }
    thetaValues_ *= constant::mathematical::pi / 180.0;
    
    if (profileDict_.found("rAxis"))
    {
        rAxis_ = vector(profileDict_.lookup("rAxis"));
    }
    else if (db().foundObject<globalOptions>("globalOptions"))
    {
        rAxis_ = db().lookupObject<globalOptions>("globalOptions").radialDirection();
    }

    calculateTheta();
 }
 
 
timeDependentAxialAzimuthalProfileFvPatchScalarField::
timeDependentAxialAzimuthalProfileFvPatchScalarField
(
    const timeDependentAxialAzimuthalProfileFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    timeDependentAxialProfileFvPatchScalarField(ptf, p, iF, mapper),
    thetaValues_(ptf.thetaValues_),
    thetaData_(ptf.thetaData_),
    thetaMethod_(interpolateTableBase::LINEAR),
    thetaTable_(timeValues_, thetaData_, thetaMethod_),
    rAxis_(ptf.rAxis_),
    theta_()
{
    calculateTheta();
}
 
 
timeDependentAxialAzimuthalProfileFvPatchScalarField::
timeDependentAxialAzimuthalProfileFvPatchScalarField
(
    const timeDependentAxialAzimuthalProfileFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    timeDependentAxialProfileFvPatchScalarField(ptf, iF),
    thetaValues_(ptf.thetaValues_),
    thetaData_(ptf.thetaData_),
    thetaMethod_(interpolateTableBase::LINEAR),
    thetaTable_(timeValues_, thetaData_, thetaMethod_),
    rAxis_(ptf.rAxis_),
    theta_()
{
    calculateTheta();
}
 
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 

void Foam::timeDependentAxialAzimuthalProfileFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Value of current time
    const scalar userTime = db().time().timeOutputValue();

    // Create scalar fields for the profile values at time t
    scalarField zData(table_(userTime));
    scalarField thetaData(thetaTable_(userTime));

    // Create interpolation table
    scalarInterpolateTable f(zValues_, zData, zMethod_);
    scalarInterpolateTable g(thetaValues_, thetaData, thetaMethod_);

    // Loop over all the face centers of the boundary patch and update
    //  the imposed patch temperature field
    forAll(z_, i)
    {
        T_[i] = f(z_[i]) + g(theta_[i]);
    }

    fixedTemperatureFvPatchScalarField::updateCoeffs();
}
 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    timeDependentAxialAzimuthalProfileFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
