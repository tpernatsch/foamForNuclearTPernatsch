/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "fixedTorqueFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

    
void fixedTorqueFvPatchVectorField::check()
{
    // Check that pressureList and tractionList from parent class are not set
    // since these will overwrite the torque BC
    if (tractionList_.valid() || pressureList_.valid())
    {
        FatalErrorInFunction() 
            << "pressureList and tractionList should not be provided." << nl
            << abort(FatalError);
    }
    
    // Check that axis can be properly normalised
    axis_ /= mag(axis_) + VSMALL;
    
    if (abs(1-mag(axis_)) > SMALL)
    {
        FatalErrorInFunction() << "Invalid axis specified " << axis_ << nl
            << abort(FatalError);
    }
}

    
void fixedTorqueFvPatchVectorField::calculateRotationCentre()
{
    // Calculate centroid
    const fvPatch& p = this->patch();
    
    const scalarField& magSf = p.magSf();
    const vectorField& Cf = p.Cf();
    
    C0_ = vector::zero;
    scalar sumArea = 0;
    
    forAll(Cf, faceI)
    {
        sumArea += magSf[faceI];
        C0_ += Cf[faceI]*magSf[faceI];
    }
    
    C0_ /= sumArea + VSMALL;
    
    if (debug)
    {
        Info<< "Rotation Centre: " << C0_ << endl;
    }
}
    
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedTorqueFvPatchVectorField::fixedTorqueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(p, iF),
    torque_(0),
    axis_(0,0,1),
    rotation_(axis_, constant::mathematical::piByTwo),
    C0_(vector::zero)
{
    calculateRotationCentre();
}


fixedTorqueFvPatchVectorField::fixedTorqueFvPatchVectorField
(
    const fixedTorqueFvPatchVectorField& pf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    tractionDisplacementFvPatchVectorField(pf, p, iF, mapper),
    torque_(pf.torque_),
    axis_(pf.axis_),
    rotation_(axis_, constant::mathematical::piByTwo),
    C0_(pf.C0_)
{}


fixedTorqueFvPatchVectorField::fixedTorqueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    tractionDisplacementFvPatchVectorField(p, iF, dict),
    torque_(readScalar(dict.lookup("torque"))),
    axis_(dict.lookup("axis")),
    rotation_(axis_, constant::mathematical::piByTwo),
    C0_(vector::zero)
{
    check();
    
    if (dict.found("centre"))
    {
        C0_ = vector(dict.lookup("centre"));
    }
    else
    {
        calculateRotationCentre();
    }
}


fixedTorqueFvPatchVectorField::fixedTorqueFvPatchVectorField
(
    const fixedTorqueFvPatchVectorField& pf,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(pf, iF),
    torque_(pf.torque_),
    axis_(pf.axis_),
    rotation_(axis_, constant::mathematical::piByTwo),
    C0_(pf.C0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedTorqueFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    tractionDisplacementFvPatchVectorField::autoMap(m);
}


void fixedTorqueFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    tractionDisplacementFvPatchVectorField::rmap(ptf, addr);
}


void fixedTorqueFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const vectorField& Sf = patch().Sf();
    const scalarField& magSf = patch().magSf();
    const vectorField& Cf = patch().Cf();
    
    scalar sumTorque = 0;
    
    forAll(Sf, faceI)
    {
        // Distance vector (radius)
        vector r = Cf[faceI] - C0_;
        r -= (r & axis_)*axis_;
        vector nr = r / mag(r);
        
        sumTorque += magSf[faceI]*mag(r);
        traction_[faceI] = rotation_.transform(nr);
    }
    
    reduce(sumTorque, sumOp<scalar>());
    
    pressure_ = 0;
    traction_ *= torque_ / (sumTorque + VSMALL);
    
    tractionDisplacementFvPatchVectorField::updateCoeffs();
}


void fixedTorqueFvPatchVectorField::write(Ostream& os) const
{
    tractionDisplacementFvPatchVectorField::write(os);

    os.writeKeyword("torque") << torque_ << token::END_STATEMENT << nl;
    os.writeKeyword("centre") << C0_ << token::END_STATEMENT << nl;
    os.writeKeyword("axis") << axis_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedTorqueFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
