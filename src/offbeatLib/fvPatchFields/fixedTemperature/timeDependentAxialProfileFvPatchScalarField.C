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

#include "timeDependentAxialProfileFvPatchScalarField.H"
#include "scalarFieldFieldINew.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#ifdef OPENFOAMFOUNDATION
#include "helperFuns.H"
#endif

namespace Foam
{

//  * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void timeDependentAxialProfileFvPatchScalarField::calculateAxialLocations()
{
    axis_ /= mag(axis_);
    z_ = referencePatch_.Cf() & axis_;

#ifdef OPENFOAMFOUNDATION
    if (debug > 1)
    {
        helperFuns::writePatchFieldVTK
        (
            patch().patch(), 
            z_,
            patch().name(),
            "z",
            patch().boundaryMesh().mesh().time().timeName()
        );
    }
#endif
}

//  * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
timeDependentAxialProfileFvPatchScalarField::
timeDependentAxialProfileFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedTemperatureFvPatchScalarField(p, iF),
    profileDict_(),
	referencePatch_(
		p.boundaryMesh().mesh().foundObject<fvMesh>("referenceMesh")
		?
		p.boundaryMesh().mesh().lookupObject<fvMesh>("referenceMesh").boundary()[patch().index()]
	  :
		p
	),
    zValues_(),
    timeValues_(),
    profileData_(),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    tMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    table_(timeValues_, profileData_, tMethod_),
	axis_(0,0,1),
    z_()
{
    if (db().foundObject<globalOptions>("globalOptions"))
    {
        axis_ = db().lookupObject<globalOptions>("globalOptions").pinDirection();
    }
    
    calculateAxialLocations();
}

timeDependentAxialProfileFvPatchScalarField::
timeDependentAxialProfileFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fixedTemperatureFvPatchScalarField(p, iF, dict, valueRequired),
	profileDict_(dict.subDict("timeDependentProfile")),
	referencePatch_(
		p.boundaryMesh().mesh().foundObject<fvMesh>("referenceMesh")
		?
		p.boundaryMesh().mesh().lookupObject<fvMesh>("referenceMesh").boundary()[patch().index()]
	  :
		p
	),
    zValues_(profileDict_.lookup("axialLocations")),
    timeValues_(profileDict_.lookup("timePoints")),
    profileData_
    (
        PtrList<scalarField>
        (
            profileDict_.lookup("axialData"), 
            scalarFieldFieldINew()
        )
    ),
    zMethod_(interpolateTableBase::interpolationMethodNames_[
        profileDict_.lookupOrDefault<word>("axialInterpolationMethod", "linear")
    ]),
    tMethod_(interpolateTableBase::interpolationMethodNames_[
        profileDict_.lookupOrDefault<word>("timeInterpolationMethod", "linear")
    ]),
    table_(timeValues_, profileData_, tMethod_),
	axis_(0,0,1),
    z_()
{
    // Set main pin direction
    if (profileDict_.found("zAxis"))
    {
        axis_ = vector(profileDict_.lookup("zAxis"));
    }
    else if (db().foundObject<globalOptions>("globalOptions"))
    {
        axis_ = db().lookupObject<globalOptions>("globalOptions").pinDirection();
    }

    calculateAxialLocations();
 }
 
 
timeDependentAxialProfileFvPatchScalarField::
timeDependentAxialProfileFvPatchScalarField
(
    const timeDependentAxialProfileFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedTemperatureFvPatchScalarField(ptf, p, iF, mapper),
    profileDict_(ptf.profileDict_),
	referencePatch_(
		p.boundaryMesh().mesh().foundObject<fvMesh>("referenceMesh")
		?
		p.boundaryMesh().mesh().lookupObject<fvMesh>("referenceMesh").boundary()[patch().index()]
	  :
		p
	),
    zValues_(ptf.zValues_),
    timeValues_(ptf.timeValues_),
    profileData_(ptf.profileData_),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    tMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    table_(timeValues_, profileData_, tMethod_),
	axis_(ptf.axis_),
    z_()
{
    calculateAxialLocations();
}
 
 
timeDependentAxialProfileFvPatchScalarField::
timeDependentAxialProfileFvPatchScalarField
(
    const timeDependentAxialProfileFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedTemperatureFvPatchScalarField(ptf, iF),
    profileDict_(ptf.profileDict_),
	referencePatch_(ptf.referencePatch_),
    zValues_(ptf.zValues_),
    timeValues_(ptf.timeValues_),
    profileData_(ptf.profileData_),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    tMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    table_(timeValues_, profileData_, tMethod_),
	axis_(ptf.axis_),
    z_()
{
    calculateAxialLocations();
}
 
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 

void Foam::timeDependentAxialProfileFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Value of current time
    const scalar userTime = db().time().timeOutputValue();

    // Create a scalar field for the profile values at time t
    scalarField zData(table_(userTime));

    // Create interpolation table     
    scalarInterpolateTable zTable(zValues_, zData, zMethod_);

    // Loop over all the face centers of the boundary patch and update
    //  the imposed patch temperature field
    forAll(z_, i)
    {
        T_[i] = zTable(z_[i]);
    }

    fixedTemperatureFvPatchScalarField::updateSurfaceTemperature();

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::timeDependentAxialProfileFvPatchScalarField::write(Ostream& os) const
{
    fixedTemperatureFvPatchScalarField::write(os);
#ifdef OPENFOAMFOUNDATION    
    writeEntry(os, "timeDependentProfile", profileDict_);
#elif OPENFOAMESI
    os.writeEntry("timeDependentProfile", profileDict_);
#endif    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    timeDependentAxialProfileFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam