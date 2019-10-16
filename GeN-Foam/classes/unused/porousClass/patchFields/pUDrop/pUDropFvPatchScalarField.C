/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "pUDropFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fieldTypes.H"
#include "rhoThermo.H"
#include "volFields.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pUDropFvPatchScalarField::pUDropFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<Foam::scalar, volMesh>& iF
)
:
    jumpCyclicFvPatchField<scalar>(p, iF),
    K_(Zero),
    UName_("U"),
    uRelax_(false),
    pInfo_(false)
{}


Foam::pUDropFvPatchScalarField::pUDropFvPatchScalarField
(
    const pUDropFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<Foam::scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    jumpCyclicFvPatchField<scalar>(ptf, p, iF, mapper),
    K_(ptf.K_),
    UName_(ptf.UName_),
    uRelax_(ptf.uRelax_),
    pInfo_(ptf.pInfo_)
{}


Foam::pUDropFvPatchScalarField::pUDropFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<Foam::scalar, volMesh>& iF,
    const dictionary& dict
)
:
    jumpCyclicFvPatchField<scalar>(p, iF),
    K_(dict.lookupOrDefault<scalar>("K", Zero)),
    UName_(dict.lookupOrDefault<word>("UName", "U")),
    uRelax_(dict.lookupOrDefault<bool>("uRelax", false)),
    pInfo_(dict.lookupOrDefault<bool>("pInfo", false))
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
        this->evaluate(Pstream::commsTypes::blocking);
    }
}


Foam::pUDropFvPatchScalarField::pUDropFvPatchScalarField
(
    const pUDropFvPatchScalarField& ptf
)
:
    jumpCyclicFvPatchField<scalar>(ptf),
    K_(ptf.K_),
    UName_(ptf.UName_),
    uRelax_(ptf.uRelax_)
{}


Foam::pUDropFvPatchScalarField::pUDropFvPatchScalarField
(
    const pUDropFvPatchScalarField& ptf,
    const DimensionedField<Foam::scalar, volMesh>& iF
)
:
    jumpCyclicFvPatchField<scalar>(ptf, iF),
    K_(ptf.K_),
    UName_(ptf.UName_),
    uRelax_(ptf.uRelax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::scalar>> Foam::pUDropFvPatchScalarField::jump() const
{
    //- Patch label to access fields at this boundary
    label patchi = this->patch().index();

    //- Access master data if needed
    scalar K(K_);
    word UName(UName_);
    bool uRelax(uRelax_);
    vectorField masterNf(this->patch().nf());
    if (!this->cyclicPatch().owner())
    {
        const fvPatch& nbrPatch =
        refCast<const fvPatch>
        (
            this->cyclicPatch().neighbFvPatch()
        );
        const pUDropFvPatchScalarField& nbrField = 
        refCast<const pUDropFvPatchScalarField>
        (
            nbrPatch.template lookupPatchField<volScalarField, scalar>(this->internalField().name())
        );
        K = nbrField.K_;
        UName = nbrField.UName_;
        uRelax = nbrField.uRelax_;
        masterNf = nbrPatch.nf();
    }

    //- Access density
    const rhoThermo& thermo = this->db().template lookupObject<rhoThermo>("thermophysicalProperties");
    scalarField rho(thermo.rho()().boundaryField()[patchi]);

    //- Access velocity fields at the patch
    volVectorField U(this->db().objectRegistry::lookupObject<volVectorField>(UName));
    vectorField Up(U.boundaryField()[patchi]);

    //- Project velocity onto patch normal
    scalarField Upn = Up & this->patch().nf();
    
    //- Compute pressure drop
    scalarField Dp = - K*rho*0.5*pow(Upn,2) * (masterNf & (Up/max(mag(Up), SMALL)));

    //- Under-relax if needed
    if (uRelax)
    {
        vectorField UpOT(U.oldTime().boundaryField()[patchi]);
        scalarField magUpn = mag(Upn);
        scalarField magUpnOT = mag(UpOT & this->patch().nf());
        scalar rf = min(min((magUpnOT)/max(magUpn, SMALL)), 1);
        Dp *= pow(rf, 2);
    }

    //- Infos
    if (pInfo_)
    {
        scalar Vdot(Zero);
        forAll(Upn, i)
        {
            Vdot += Upn[i]*this->patch().magSf()[i];
        }
        Info << "Patch " << this->patch().name() << ", avg(Dp) = " << gAverage(Dp) << " Pa, Vdot = " << Vdot << " m3/s "<< endl;
    }

    tmp<scalarField> tDp(new scalarField(Dp));
    return tDp;
}


void Foam::pUDropFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    jumpCyclicFvPatchField<scalar>::autoMap(m);
}


void Foam::pUDropFvPatchScalarField::rmap
(
    const fvPatchField<Foam::scalar>& ptf,
    const labelList& addr
)
{
    jumpCyclicFvPatchField<scalar>::rmap(ptf, addr);
}


void Foam::pUDropFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("patchType") << this->interfaceFieldType()
        << token::END_STATEMENT << nl;

    if (this->cyclicPatch().owner())
    {
        os.writeKeyword("K")<< K_
            << token::END_STATEMENT << nl;
        os.writeKeyword("UName")<< UName_
            << token::END_STATEMENT << nl;
        os.writeKeyword("uRelax")<< uRelax_
            << token::END_STATEMENT << nl;
        os.writeKeyword("pInfo")<< pInfo_
            << token::END_STATEMENT << nl;
    }

    this->writeEntry("value", os);
}


// ************************************************************************* //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        pUDropFvPatchScalarField
    );
}