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

#include "topCladRingPressureFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

#include "gapGasModel.H"
#include "globalOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

topCladRingPressureFvPatchVectorField::
topCladRingPressureFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    coolantPressureFvPatchVectorField(p, iF),
    topCapRadius_(0.0),
    topCapPlenumRadius_(0.0),
    fuelTopPatchID_(-1),
    displacementName_("D")
{}


topCladRingPressureFvPatchVectorField::
topCladRingPressureFvPatchVectorField
(
    const topCladRingPressureFvPatchVectorField& tvtdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coolantPressureFvPatchVectorField(tvtdpvf, p, iF, mapper),
    topCapRadius_(tvtdpvf.topCapRadius_),
    topCapPlenumRadius_(tvtdpvf.topCapPlenumRadius_),
    fuelTopPatchID_(tvtdpvf.fuelTopPatchID_),
    displacementName_(tvtdpvf.displacementName_)
{}


topCladRingPressureFvPatchVectorField::
topCladRingPressureFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    coolantPressureFvPatchVectorField(p, iF, dict),
    topCapRadius_(readScalar(dict.lookup("topCapOuterRadius"))),
    topCapPlenumRadius_(readScalar(dict.lookup("topCapInnerRadius"))),
    fuelTopPatchID_(-1),
    displacementName_(dict.lookupOrDefault<word>("displacementName", "D"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void topCladRingPressureFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    coolantPressureFvPatchVectorField::autoMap(m);
}


void topCladRingPressureFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    coolantPressureFvPatchVectorField::rmap(ptf, addr);
}


void topCladRingPressureFvPatchVectorField::updateTraction()
{
    // Set fuelTopPatchID_
    if(fuelTopPatchID_ == -1)
    {
        const volVectorField& vfD =
            this->db().lookupObject<volVectorField>(internalField().name());

        int nFuelTopPatch(0);

        forAll(vfD.boundaryField(), patchID)
        {
            if(isA<plenumSpringBase>(vfD.boundaryField()[patchID]))
            {
                fuelTopPatchID_ = patchID;
                nFuelTopPatch += 1;
            }
        }

        if(nFuelTopPatch > 1)
        {
            FatalErrorInFunction()
            << "More than one plenumSpring type patch was found." << nl
            << "topCladRingPressure patchField requires one and "
            << "only one plenum spring patch-type in the model." 
            << abort(FatalError);
        }

        if(nFuelTopPatch == 0)
        {
            FatalErrorInFunction()
            << "No plenumSpring type patch was found." << nl
            << "topCladRingPressure patchField requires one and "
            << "only one plenum spring patch-type in the model." 
            << abort(FatalError);
        }
    }

    // Pi
    scalar pi = Foam::constant::mathematical::pi;

    const globalOptions& globalOpt
    (db().lookupObject<globalOptions>("globalOptions"));

    scalar angularFraction = globalOpt.angularFraction();

    // Include the gap gas pressure (only when the cap is not explicitly modeled)
    const gapGasModel& gapGas
    = db().lookupObject<gapGasModel>("gapGas");

    // Reference to total displacement field
    const volVectorField& D =
        db().lookupObject<volVectorField>(internalField().name());

    // Request spring elongation and spring modulus from fuel top patch
    const plenumSpringBase& fuelTopD(
        refCast<const plenumSpringBase>(D.boundaryField()[fuelTopPatchID_]));

    scalar Dtot(fuelTopD.springElongation());
    scalar springModulus(fuelTopD.springModulus());

    // Top clad ring area
    scalar Atot = gSum(patch().magSf())/angularFraction;

    // Top cap area fluid side. 
    scalar topCapArea = pi*pow(topCapRadius_, 2.0);

    // Top cap area plenum side. 
    scalar topCapPlenumArea = pi*pow(topCapPlenumRadius_, 2.0);
      
    // Set pressures

    // Inner gap gas pressure (negative because it acts in the opposite direction
    // on the cladding top ring)
    scalarField gapGasPressure
    (patch().size(), -gapGas.p()*topCapPlenumArea/Atot);
      
    // Spring pressure (negative because it acts in the opposite direction
    // on the cladding top ring)
    scalarField springPressure
    (patch().size(), -springModulus*Dtot/Atot);

    // Fluid pressure
    if(coolantPressureList_.valid())
    {
        coolantPressure_ = 
        coolantPressureList_->value
        (
            this->db().time().timeToUserTime
            (
                this->db().time().value()
            )
        );
    }

    // NOTE: cannot multiply directly coolantPressure_ --> create new scalarField
    scalarField coolantPressure(coolantPressure_*topCapArea/Atot);      

    // Update total pressure
    pressure_ = (gapGasPressure + springPressure + coolantPressure);
}


void topCladRingPressureFvPatchVectorField::write(Ostream& os) const
{
    coolantPressureFvPatchVectorField::write(os);

#ifdef OPENFOAMFOUNDATION  
    writeEntry<scalar>(os, "topCapOuterRadius", topCapRadius_);
    writeEntry<scalar>(os, "topCapInnerRadius", topCapPlenumRadius_);

    writeEntryIfDifferent<word>(os, "displacementName", "D", displacementName_);
#elif OPENFOAMESI
    os.writeEntry<scalar>("topCapOuterRadius", topCapRadius_);
    os.writeEntry<scalar>("topCapInnerRadius", topCapPlenumRadius_);

    os.writeEntryIfDifferent<word>("displacementName", "D", displacementName_);
#endif

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    topCladRingPressureFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
