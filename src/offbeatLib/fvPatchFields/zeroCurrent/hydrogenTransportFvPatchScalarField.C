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

#include "hydrogenTransportFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "hydrogenTransportSolver.H"
#ifdef OPENFOAMFOUNDATION
#include "helperFuns.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hydrogenTransportFvPatchScalarField::
hydrogenTransportFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchScalarField(p, iF),
    ingressRate_(p.size(), 0),
    phi_(p.size(), 0),
    deltaCoeffsD_(p.size(), 0)
{}


hydrogenTransportFvPatchScalarField::
hydrogenTransportFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchScalarField(p, iF),
    ingressRate_("ingressRate", dict, p.size()),
    phi_(p.size(), 0),
    deltaCoeffsD_(p.size(), 0)
{}


hydrogenTransportFvPatchScalarField::
hydrogenTransportFvPatchScalarField
(
    const hydrogenTransportFvPatchScalarField& tdpvf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchScalarField(tdpvf, p, iF, mapper),
    ingressRate_(tdpvf.ingressRate_),
    phi_(tdpvf.phi_),
    deltaCoeffsD_(tdpvf.deltaCoeffsD_)
{}    


hydrogenTransportFvPatchScalarField::
hydrogenTransportFvPatchScalarField
(
    const hydrogenTransportFvPatchScalarField& tdpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:  
    fvPatchScalarField(tdpvf, iF),
    ingressRate_(tdpvf.ingressRate_),
    phi_(tdpvf.phi_),
    deltaCoeffsD_(tdpvf.deltaCoeffsD_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void hydrogenTransportFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchScalarField::autoMap(m);

#ifdef OPENFOAMFOUNDATION
    m(ingressRate_, ingressRate_);
    m(phi_, phi_);
    m(deltaCoeffsD_, deltaCoeffsD_);
#elif OPENFOAMESI
    ingressRate_.autoMap(m);
    phi_.autoMap(m);
    deltaCoeffsD_.autoMap(m);
#endif
    
}


void hydrogenTransportFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fvPatchScalarField::rmap(ptf, addr);

    const hydrogenTransportFvPatchScalarField& fgptf =
        refCast<const hydrogenTransportFvPatchScalarField>(ptf);

    ingressRate_.rmap(fgptf.ingressRate_, addr);
    phi_.rmap(fgptf.phi_, addr);
    deltaCoeffsD_.rmap(fgptf.deltaCoeffsD_, addr);
}


void hydrogenTransportFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const hydrogenTransportSolver& model = 
        this->patch().boundaryMesh().mesh().lookupObject<hydrogenTransportSolver>
        (
            hydrogenTransportSolver::typeName
        );
    
    const word& diffusionCoeffName = model.diffH().name();
    const word& TName = model.T().name();
    //scalar Q = model.Q().value();
        
    //- diffusion coefficient patch
    const fvPatchScalarField& diffHf =
        patch().lookupPatchField<volScalarField, scalar>(diffusionCoeffName);

    // Reference to temperature gradient patch field
    const fvPatchScalarField& Tf = 
        patch().lookupPatchField<volScalarField, scalar>(TName);
    scalarField snGradTf = Tf.snGrad();

    // Build Soret term
     //scalar R = constant::physicoChemical::R.value();
     //phi_ = diffHf*Q / (R*pow(Tf, 2))*snGradTf;

    //const transportSolver& stressDrivenDiff = this->patch().boundaryMesh().mesh().lookupObject<transportSolver>(HNGD::typeName);
    //bool StressDrivenDiffusion_ = stressDrivenDiff.lookupOrDefault<bool>("stressDrivenDiffusion", false);
    /*
    if(stressDrivenDiffusion)
    {
        phi_ -= diffHf * V_H / (R*Tf)* snGradsigmaHydf;
    }
    */

    // Reference to flux patch field
    //Flux (phi) from Soret term and Stress diffusion (if present)
    phi_ = patch().lookupPatchField<surfaceScalarField, scalar>("phiH") / patch().magSf();

    deltaCoeffsD_ = diffHf * patch().deltaCoeffs();
    
    fvPatchScalarField::updateCoeffs();
}


void hydrogenTransportFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

#ifdef OPENFOAMFOUNDATION
    if (debug > 1)
    {
        helperFuns::writePatchFieldVTK
        (
            patch().patch(), 
            deltaCoeffsD_,
            patch().name(),
            "deltaCoeffsD",
            patch().boundaryMesh().mesh().time().timeName()
        );
        
        helperFuns::writePatchFieldVTK
        (
            patch().patch(), 
            phi_,
            patch().name(),
            "phi",
            patch().boundaryMesh().mesh().time().timeName()
        );
    }
#endif
    
    fvPatchScalarField::operator=(
        (deltaCoeffsD_*patchInternalField() + ingressRate_)
      / (deltaCoeffsD_ + phi_)
    );
    
    fvPatchScalarField::evaluate();
}


tmp<scalarField>
hydrogenTransportFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return deltaCoeffsD_ / (deltaCoeffsD_ + phi_);
    //return tmp<scalarField>(new scalarField(size(), 0.0));
}


tmp<scalarField>
hydrogenTransportFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return ingressRate_ / (deltaCoeffsD_ + phi_);
    //return tmp<scalarField>(new scalarField(size(), 0.0));
}


tmp<scalarField>
hydrogenTransportFvPatchScalarField::gradientInternalCoeffs() const
{
    const scalarField& deltaCoeffs = patch().deltaCoeffs();

    return -phi_*deltaCoeffs / (deltaCoeffsD_ + phi_);
    //return tmp<scalarField>(new scalarField(size(), 0.0));
}


tmp<scalarField>
hydrogenTransportFvPatchScalarField::gradientBoundaryCoeffs() const
{
    const scalarField& deltaCoeffs = patch().deltaCoeffs();

    return ingressRate_*deltaCoeffs / (deltaCoeffsD_ + phi_);
    //return ingressRate_*deltaCoeffs / deltaCoeffsD_;
}


void hydrogenTransportFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
#ifdef OPENFOAMFOUNDATION
    writeEntry(os, "ingressRate", ingressRate_);
    writeEntry(os, "value", *this);
#elif OPENFOAMESI    
    os.writeEntry("ingressRate", ingressRate_);
    os.writeEntry("value", *this);
#endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    hydrogenTransportFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
