/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2506                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2025 OpenCFD Ltd.          |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "convectiveFromFluidFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatch.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Time.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * //

convectiveFromFluidFvPatchScalarField::convectiveFromFluidFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    kappaName_("kappa"),
    htcName_("mappedHtc"),
    TfluidName_("mappedTFluid"),
    TFluid_
    (
        IOobject
        (
            "mappedTFluid",
            this->patch().boundaryMesh().mesh().time().constant(),
            this->patch().boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->patch().boundaryMesh().mesh(),
        dimensionedScalar("mappedTFluid", dimTemperature, 0.0)
    ),
    htc_
    (
        IOobject
        (
            "mappedHtc",
            this->patch().boundaryMesh().mesh().time().constant(),
            this->patch().boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->patch().boundaryMesh().mesh(),
        dimensionedScalar("mappedHtc", dimPower/pow(dimLength,2)/dimTemperature, 0.0)
    )
{
    // Default: Dirichlet to zero
    refValue() = 0.0;
    refGrad()  = 0.0;
    valueFraction() = 1.0;
}


convectiveFromFluidFvPatchScalarField::convectiveFromFluidFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict),
    kappaName_(dict.lookupOrDefault<word>("kappa", "k")),
    htcName_(dict.lookupOrDefault<word>("htc", "mappedHtc")),
    TfluidName_(dict.lookupOrDefault<word>("Tfluid", "mappedTFluid")),
    TFluid_
    (
        IOobject
        (
            "mappedTFluid",
            this->patch().boundaryMesh().mesh().time().constant(),
            this->patch().boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->patch().boundaryMesh().mesh(),
        dimensionedScalar("mappedTFluid", dimTemperature, 0.0)
    ),
    htc_
    (
        IOobject
        (
            "mappedHtc",
            this->patch().boundaryMesh().mesh().time().constant(),
            this->patch().boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->patch().boundaryMesh().mesh(),
        dimensionedScalar("mappedHtc", dimPower/pow(dimLength,2)/dimTemperature, 0.0)
    )
{
    if (dict.found("refValue")) refValue() = scalarField("refValue", dict, p.size());
    if (dict.found("refGradient")) refGrad() = scalarField("refGradient", dict, p.size());
    if (dict.found("valueFraction")) valueFraction() = scalarField("valueFraction", dict, p.size());
}

convectiveFromFluidFvPatchScalarField::convectiveFromFluidFvPatchScalarField
(
    const convectiveFromFluidFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),  // <- only ptf and iF
    kappaName_(ptf.kappaName_),
    htcName_(ptf.htcName_),
    TfluidName_(ptf.TfluidName_),
    TFluid_(ptf.TFluid_),
    htc_(ptf.htc_)
{}

convectiveFromFluidFvPatchScalarField::convectiveFromFluidFvPatchScalarField
(
    const convectiveFromFluidFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& m
)
:
    mixedFvPatchScalarField(ptf, p, iF, m),
    kappaName_(ptf.kappaName_),
    htcName_(ptf.htcName_),
    TfluidName_(ptf.TfluidName_),
    TFluid_(ptf.TFluid_),
    htc_(ptf.htc_)
{}


// * * * * * * * * * * * * * * * * Member Functions * * * * * * * * * //

void convectiveFromFluidFvPatchScalarField::updateCoeffs()
{
    if (updated()) return;

    const volScalarField& kappa =
        this->patch().boundaryMesh().mesh().lookupObject<volScalarField>(kappaName_);

    const fvPatch& p = patch();
    const labelList& faceCells = p.faceCells();

    forAll(p, faceI)
    {
        label cellI = faceCells[faceI];

        scalar h = htc_[cellI];
        scalar Tf = TFluid_[cellI];
        scalar k = kappa[cellI];

        // distance coefficient (deltaCoeffs = 1/d)
        scalar d = 1.0/p.deltaCoeffs()[faceI];

        refValue()[faceI] = Tf;
        refGrad()[faceI]  = 0.0;
        valueFraction()[faceI] = h / (h + k/d);
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void convectiveFromFluidFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("kappa") << kappaName_ << token::END_STATEMENT << nl;
    os.writeKeyword("htc") << htcName_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tfluid") << TfluidName_ << token::END_STATEMENT << nl;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "addToRunTimeSelectionTable.H"
namespace Foam
{
    makePatchTypeField(fvPatchScalarField, convectiveFromFluidFvPatchScalarField);
}
