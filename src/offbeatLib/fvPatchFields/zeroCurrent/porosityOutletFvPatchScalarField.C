 /*---------------------------------------------------------------------------*\
   =========                 |
   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     | Website:  https://openfoam.org
     \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
 
 #include "porosityOutletFvPatchScalarField.H"
 #include "addToRunTimeSelectionTable.H"
 #include "fvPatchFieldMapper.H"
 #include "volFields.H"
 #include "EulerDdtScheme.H"
 #include "CrankNicolsonDdtScheme.H"
 #include "backwardDdtScheme.H"
 #include "localEulerDdtScheme.H"

 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
 
 
 // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
 
 Foam::porosityOutletFvPatchScalarField::porosityOutletFvPatchScalarField
 (
     const fvPatch& p,
     const DimensionedField<scalar, volMesh>& iF
 )
 :
     mixedFvPatchField<scalar>(p, iF),
     factor_(1.0)
 {
     this->refValue() = Zero;
     this->refGrad() = Zero;
     this->valueFraction() = 0.0;
 }
 
 
 Foam::porosityOutletFvPatchScalarField::porosityOutletFvPatchScalarField
 (
     const fvPatch& p,
     const DimensionedField<scalar, volMesh>& iF,
     const dictionary& dict
 )
 :
     mixedFvPatchField<scalar>(p, iF),
     factor_(dict.lookupOrDefault("factor", 1.0))
 {
     if (dict.found("value"))
     {
         fvPatchField<scalar>::operator=
         (
             Field<scalar>("value", dict, p.size())
         );
     }
     else
     {
         fvPatchField<scalar>::operator=(this->patchInternalField());
     }
 
     this->refValue() = *this;
     this->refGrad() = Zero;
     this->valueFraction() = 0.0;

 }
 
 
 Foam::porosityOutletFvPatchScalarField::porosityOutletFvPatchScalarField
 (
     const porosityOutletFvPatchScalarField& ptf,
     const fvPatch& p,
     const DimensionedField<scalar, volMesh>& iF,
     const fvPatchFieldMapper& mapper
 )
 :
     mixedFvPatchField<scalar>(ptf, p, iF, mapper),
     factor_(ptf.factor_)
 {}
 
 
 Foam::porosityOutletFvPatchScalarField::porosityOutletFvPatchScalarField
 (
     const porosityOutletFvPatchScalarField& ptpsf,
     const DimensionedField<scalar, volMesh>& iF
 )
 :
     mixedFvPatchField<scalar>(ptpsf, iF),
     factor_(ptpsf.factor_)
 {}
 
 
 // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
 scalarField Foam::porosityOutletFvPatchScalarField::advectionSpeed() const
 {
    // - METHOD 1: lookup vCoalescence_
    if(db().foundObject<volVectorField>("vCoalescence"))
    {
        const fvPatchVectorField& coalVelP =
        this->patch().lookupPatchField<volVectorField, vector>
        (
            "vCoalescence"
        );
        
        return coalVelP & this->patch().nf();
    }
    else
    {
        return scalarField(this->size(), 0.0);
    }


    // - METHOD 2: based on deltaPorosity
    // if(db().foundObject<volScalarField>("porosity"))
    // {
    //     // const fvPatchScalarField& poro =
    //     // this->patch().lookupPatchField<volScalarField, scalar>
    //     // (
    //     //     "porosity"
    //     // );

    //     const volScalarField& porosity =
    //         db().lookupObject<volScalarField>("porosity");

    //     scalarField poro(porosity, patch().faceCells());

    //     const scalar& deltaT = db().time().deltaTValue();
        
    //     scalarField thrPoro(this->size(), 0.75);

    //     const fvMesh& mesh = patch().boundaryMesh().mesh();

    //     scalarField V(mesh.V(), patch().faceCells());

    //     scalarField vel((max(poro-thrPoro, 0.0)*V)/(patch().magSf()*deltaT));

    //     return factor_ * vel;
    // }
    // else
    // {
    //     return scalarField(this->size(), 0.0);
    // }
 }
 
 
 void Foam::porosityOutletFvPatchScalarField::updateCoeffs()
 {
    if (this->updated())
    {
        return;
    }

    const fvMesh& mesh = this->internalField().mesh();

    word ddtScheme
    (
        mesh.ddtScheme(this->internalField().name())
    );
    scalar deltaT = this->db().time().deltaTValue();

    const volScalarField& field = this->db().lookupObject<volScalarField>
    (
        this->internalField().name()
    );

    // Calculate the advection speed of the field wave
    // If the wave is incoming set the speed to 0.
    const scalarField w(Foam::max(advectionSpeed(), scalar(0)));

    // Calculate the field wave coefficient alpha (See notes)
    const scalarField alpha(w*deltaT*this->patch().deltaCoeffs());

    label patchi = this->patch().index();
 
    
    if
    (
        ddtScheme == fv::EulerDdtScheme<scalar>::typeName
    || ddtScheme == fv::CrankNicolsonDdtScheme<scalar>::typeName
    )
    {
        this->refValue() = field.oldTime().boundaryField()[patchi];

        this->valueFraction() = 1.0/(1.0 + alpha);
    }
    else if (ddtScheme == fv::backwardDdtScheme<scalar>::typeName)
    {
        this->refValue() =
        (
            2.0*field.oldTime().boundaryField()[patchi]
        - 0.5*field.oldTime().oldTime().boundaryField()[patchi]
        )/1.5;

        this->valueFraction() = 1.5/(1.5 + alpha);
    }
    else if
    (
        ddtScheme == fv::localEulerDdtScheme<scalar>::typeName
    )
    {
        const volScalarField& rDeltaT =
            fv::localEulerDdt::localRDeltaT(mesh);

        // Calculate the field wave coefficient alpha (See notes)
        const scalarField alpha
        (
            w*this->patch().deltaCoeffs()/rDeltaT.boundaryField()[patchi]
        );

        this->refValue() = field.oldTime().boundaryField()[patchi];

        this->valueFraction() = 1.0/(1.0 + alpha);
    }
    else
    {
        FatalErrorInFunction
            << ddtScheme
            << "\n    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }
 
    mixedFvPatchField<scalar>::updateCoeffs();
 }
 
 
 void Foam::porosityOutletFvPatchScalarField::write(Ostream& os) const
 {
    fvPatchField<scalar>::write(os);

#ifdef OPENFOAMFOUNDATION
    writeEntry(os, "velocity", advectionSpeed());
    writeEntry(os, "value", *this);
#elif OPENFOAMESI    
    advectionSpeed().writeEntry("velocity", os);
    this->writeEntry("value", os);
#endif    
 }
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    porosityOutletFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}
 // ************************************************************************* //