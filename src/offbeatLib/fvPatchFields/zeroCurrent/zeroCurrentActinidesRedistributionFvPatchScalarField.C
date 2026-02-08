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

#include "zeroCurrentActinidesRedistributionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

zeroCurrentActinidesRedistributionFvPatchScalarField::
zeroCurrentActinidesRedistributionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    R_(8.31446261815324),
    Q_(-146.5e3),
    poreDiameter_(0.008e-3),
    poreThickness_(80e-3),
    A_(0.35)
{
}

zeroCurrentActinidesRedistributionFvPatchScalarField::
zeroCurrentActinidesRedistributionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    R_(8.31446261815324),
    Q_(-146.5e3),
    poreDiameter_(0.008e-3),
    poreThickness_(80e-3),
    A_(0.35)
{

    fvPatchScalarField::operator=(patchInternalField());
    gradient() = 0;
}

zeroCurrentActinidesRedistributionFvPatchScalarField::
zeroCurrentActinidesRedistributionFvPatchScalarField
(
    const zeroCurrentActinidesRedistributionFvPatchScalarField& tdpvf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(tdpvf, p, iF, mapper),
    R_(8.31446261815324),
    Q_(-146.5e3),
    poreDiameter_(0.008e-3),
    poreThickness_(80e-3),
    A_(0.35)
{
}    


zeroCurrentActinidesRedistributionFvPatchScalarField::
zeroCurrentActinidesRedistributionFvPatchScalarField
(
    const zeroCurrentActinidesRedistributionFvPatchScalarField& tdpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:  
    fixedGradientFvPatchScalarField(tdpvf, iF),
    R_(tdpvf.R_),
    Q_(tdpvf.Q_),
    poreDiameter_(tdpvf.poreDiameter_),
    poreThickness_(tdpvf.poreThickness_),
    A_(tdpvf.A_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void zeroCurrentActinidesRedistributionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
}


void zeroCurrentActinidesRedistributionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);
}


void zeroCurrentActinidesRedistributionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    // Patch normal fields 
    const vectorField n(this->patch().nf());

    // Reference to temperature patch field
    const fvPatchScalarField& Tf =
        this->patch().lookupPatchField<volScalarField, scalar>("T");

    // Reference to temperature gradient patch field
    const fvPatchVectorField& gradTf =
        this->patch().lookupPatchField<volVectorField, vector>("gradT");

    // // Reference to porosity patch field
    // const fvPatchScalarField& pf =
    //     this->patch().lookupPatchField<volScalarField, scalar>("porosity");

    // // Reference to pore velocity
    // const fvPatchVectorField& vf =
    //     patch().lookupPatchField<volVectorField, vector>("poreVelocity");

    // Get name of the field to which this BC is applied (Pu or Am?)
    const word& fieldName = internalField().name();

    // Build diffusion coeff name as diff+fieldName
    const word diffName = "diff"+fieldName.substr(0,2);

    // Reference to Diffusion coeff
    const fvPatchField<scalar>& diff 
        = patch().lookupPatchField<volScalarField, vector>(diffName);

    // Build Soret term
    const vectorField soretTerm = Q_/(R_*pow(Tf,2)) * gradTf;

    // // Build pore migration term
    // vectorField pmTerm(soretTerm.size(), vector::zero);
    
    // Init gradient 
    scalarField g = ((*this * (*this-1) * soretTerm ) & n); 
    
    // if( max( mag(vf) ) > 1e-10 )
    // {
    //     pmTerm = pf/poreDiameter_ * A_*poreThickness_* gradTf
    //             * exp(-diff/poreThickness_/max(mag(vf),VSMALL));

    //     g += ((*this * pmTerm ) & n);
    // }
    
    // Set the gradient 
    gradient() = g; 

    fixedGradientFvPatchScalarField::updateCoeffs();
}

void zeroCurrentActinidesRedistributionFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
#ifdef OPENFOAMFOUNDATION
    writeEntry(os, "value", *this);
#elif OPENFOAMESI    
    this->writeEntry("value", os);
#endif    

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    zeroCurrentActinidesRedistributionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
