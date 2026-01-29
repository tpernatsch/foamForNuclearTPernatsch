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

#include "fixedTemperatureFvPatchScalarField.H"
#include "scalarFieldFieldINew.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{
 
//* * * * * * * * * * * * * Protected member functions  * * * * * * * * * * *//
 
void Foam::fixedTemperatureFvPatchScalarField::updateSurfaceTemperature()
{
    if
    (
        patch().boundaryMesh().mesh().foundObject<surfaceScalarField>
        ("oxideThickness")
    )
    {
        //- Current patch temperature 
        scalarField Tp(*this);

        //- Oxide thicnkess
        const scalarField& S
        = patch().lookupPatchField<surfaceScalarField, scalar>("oxideThickness");

        //- Internal field
        const scalarField Ti(patchInternalField());
        
        //- Material conductivity
        const scalarField& k
        = patch().lookupPatchField<volScalarField, scalar>(kappaName_);

        //- Compute oxide conductivity (Matpro)
        const scalarField kox = 0.835 + 1.81e-4 * Tp; 

        //- Inverse cell to patch face distances
        const scalarField d = patch().deltaCoeffs();  

        //- Thermal resistances
        scalarField R(1/d/k);
        scalarField Rox(S/kox);
        scalarField Reff(R + Rox);

        //- Loop over all the faces of the boundary patch
        forAll(patch(), i)
        {    
            //- Compute Fourier heat flux
            const scalar q = -(T_[i] - Ti[i])/Reff[i];

            //- Compute outer cladding or metal/oxide interface temperature
            Tp[i] = (T_[i] + q*S[i]/kox[i]);
        }
        
        fvPatchField<scalar>::operator==(Tp);
    }
    else
    {
        
        fvPatchField<scalar>::operator==(T_);
    }
}


//  * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //
 
fixedTemperatureFvPatchScalarField::fixedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    T_(*this),
    kappaName_("k")
{}

fixedTemperatureFvPatchScalarField::fixedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict, false),
    T_(*this),
    kappaName_(dict.lookupOrDefault<word>("kappa", "k"))
{
    if(valueRequired)
    {
        T_ = scalarField("value", dict, p.size());

        // If outer oxide temperature is specified, then the fixed T_ field is 
        // equal to materialOxideInterfaceTemperature
        if(dict.found("materialOxideInterfaceTemperature"))
        {
            fvPatchField<scalar>::operator==
            (
                scalarField
                (
                    "materialOxideInterfaceTemperature",
                    dict,
                    p.size()
                )
            );
        }
        // Otherwise initial value of T_ is equal to "value"
        else
        {
            fvPatchField<scalar>::operator==(T_);
        }
    }
}

 
fixedTemperatureFvPatchScalarField::fixedTemperatureFvPatchScalarField
(
    const fixedTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    T_(mapper(ptf.T_)),
    kappaName_(ptf.kappaName_)
{}
 
 
fixedTemperatureFvPatchScalarField::
fixedTemperatureFvPatchScalarField
(
    const fixedTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    T_(ptf.T_),
    kappaName_(ptf.kappaName_)
{}
 

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::fixedTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<scalar>::autoMap(m);
    
#ifdef OPENFOAMFOUNDATION
    m(T_, T_);
#elif OPENFOAMESI
    T_.autoMap(m);
#endif
}


void Foam::fixedTemperatureFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<scalar>::rmap(ptf, addr);

    const fixedTemperatureFvPatchScalarField& mptf =
     refCast<const fixedTemperatureFvPatchScalarField>(ptf);

    T_.rmap(mptf.T_, addr);
}


void Foam::fixedTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    updateSurfaceTemperature();

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::fixedTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    #ifdef OPENFOAMFOUNDATION
    writeEntryIfDifferent<scalarField>(os, 
        "materialOxideInterfaceTemperature", T_, *this);

    writeEntry(os, "value", T_);
    
    writeEntryIfDifferent<word>(os, "kappa", "k", kappaName_);
    
    #elif OPENFOAMESI
    os.writeEntryIfDifferent<scalarField>(
            "materialOxideInterfaceTemperature", T_, *this);

    T_.writeEntry("value", os);
    
    os.writeEntryIfDifferent<word>("kappa", "k", kappaName_);
    #endif    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    fixedTemperatureFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam