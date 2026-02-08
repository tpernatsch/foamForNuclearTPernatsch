/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "axialProfileHeatFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "scalarFieldFieldINew.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//  * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void axialProfileHeatFluxFvPatchScalarField::calculateAxialLocations()
{
    axis_ /= mag(axis_);
    z_ = referencePatch_.Cf() & axis_;
    
    //- Renormalise from 0-1
    const scalar zMin = gMin(z_);
    const scalar zMax = gMax(z_);
    z_ = (z_ - zMin) / (zMax - zMin);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

axialProfileHeatFluxFvPatchScalarField::axialProfileHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedHeatFluxFvPatchScalarField(p, iF),
    profileDict_(),
	referencePatch_(
		p.boundaryMesh().mesh().foundObject<fvMesh>("referenceMesh")
		?
		p.boundaryMesh().mesh().lookupObject<fvMesh>("referenceMesh").boundary()[patch().index()]
	  :
		p
	),
    timePoints_(),
    axialPoints_(),
    averageData_(),
    axialData_(),
    tMethod_(interpolationMethod::LINEAR),
    zMethod_(interpolationMethod::LINEAR),
    averageTable_(timePoints_, averageData_, tMethod_),
    axialTable(timePoints_, axialData_, tMethod_),
	axis_(0,0,1),
    z_()
{
    if (db().foundObject<globalOptions>("globalOptions"))
    {
        axis_ = db().lookupObject<globalOptions>("globalOptions").pinDirection();
    }
    
    calculateAxialLocations();
}


axialProfileHeatFluxFvPatchScalarField::axialProfileHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedHeatFluxFvPatchScalarField(p, iF, dict),
	profileDict_(dict.subDict("timeDependentProfile")),
	referencePatch_(
		p.boundaryMesh().mesh().foundObject<fvMesh>("referenceMesh")
		?
		p.boundaryMesh().mesh().lookupObject<fvMesh>("referenceMesh").boundary()[patch().index()]
	  :
		p
	),
    timePoints_(profileDict_.lookup("timePoints")),
    axialPoints_(profileDict_.lookup("axialLocations")),
    averageData_(profileDict_.lookup("averageData")),
    axialData_
    (
        PtrList<scalarField>
        (
            profileDict_.lookup("axialData"), 
            scalarFieldFieldINew()
        )
    ),
    tMethod_
    (
        interpolateTableBase::interpolationMethodNames_
        [
            profileDict_.lookupOrDefault<word>
            ("timeInterpolationMethod", "linear")
        ]
    ),
    zMethod_
    (
        interpolateTableBase::interpolationMethodNames_
        [
            profileDict_.lookupOrDefault<word>
            ("axialInterpolationMethod", "linear")
        ]
    ),
    averageTable_(timePoints_, averageData_, tMethod_),
    axialTable(timePoints_, axialData_, tMethod_),
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


axialProfileHeatFluxFvPatchScalarField::axialProfileHeatFluxFvPatchScalarField
(
    const axialProfileHeatFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedHeatFluxFvPatchScalarField(ptf, p, iF, mapper),
    profileDict_(ptf.profileDict_),
	referencePatch_(ptf.referencePatch_),
    timePoints_(ptf.timePoints_),
    axialPoints_(ptf.axialPoints_),
    averageData_(ptf.averageData_),
    axialData_(ptf.axialData_),
    tMethod_(ptf.tMethod_),
    zMethod_(ptf.zMethod_),
    averageTable_(timePoints_, averageData_, tMethod_),
    axialTable(timePoints_, axialData_, tMethod_),
	axis_(ptf.axis_),
    z_()
{
    calculateAxialLocations();
}


axialProfileHeatFluxFvPatchScalarField::axialProfileHeatFluxFvPatchScalarField
(
    const axialProfileHeatFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedHeatFluxFvPatchScalarField(ptf, iF),
    profileDict_(ptf.profileDict_),
	referencePatch_(ptf.referencePatch_),
    timePoints_(ptf.timePoints_),
    axialPoints_(ptf.axialPoints_),
    averageData_(ptf.averageData_),
    axialData_(ptf.axialData_),
    tMethod_(ptf.tMethod_),
    zMethod_(ptf.zMethod_),
    averageTable_(timePoints_, averageData_, tMethod_),
    axialTable(timePoints_, axialData_, tMethod_),
	axis_(ptf.axis_),
    z_()
{
    calculateAxialLocations();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void axialProfileHeatFluxFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedHeatFluxFvPatchScalarField::autoMap(m);
}


void axialProfileHeatFluxFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedHeatFluxFvPatchScalarField::rmap(ptf, addr);
}


void axialProfileHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    
    //- Interpolate in time
    const scalar time = this->db().time().value();
    const scalar userTime = this->db().time().timeToUserTime(time);
    const scalar qAvg = averageTable_(userTime);
    scalarField profile = axialTable(userTime);
    scalarInterpolateTable f(axialPoints_, profile, zMethod_);
    
    // Interpolate in axial direction
    forAll(z_, i)
    {
        q_[i] = qAvg*f(z_[i]);
    }

    fixedHeatFluxFvPatchScalarField::updateCoeffs();
}


void axialProfileHeatFluxFvPatchScalarField::write(Ostream& os) const
{
    fixedHeatFluxFvPatchScalarField::write(os);
#ifdef OPENFOAMFOUNDATION    
    writeEntry(os, "timeDependentProfile", profileDict_);
#elif OPENFOAMESI
    os.writeEntry("timeDependentProfile", profileDict_);
#endif    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, axialProfileHeatFluxFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
