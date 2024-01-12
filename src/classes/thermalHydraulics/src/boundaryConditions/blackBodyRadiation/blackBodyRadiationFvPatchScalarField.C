/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2306                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2022 OpenCFD Ltd.         |
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

#include "blackBodyRadiationFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "constants.H"
#include "fvCFD.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blackBodyRadiationFvPatchScalarField::
blackBodyRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    emissivity_(Zero),
    kappa_(Zero),
    Ta_(Zero)
{
}


Foam::blackBodyRadiationFvPatchScalarField::
blackBodyRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    emissivity_(dict.get<scalar>("emissivity")),
    kappa_(dict.get<scalar>("kappa")),
    Ta_(dict.get<scalar>("Ta"))
{
    fixedValueFvPatchScalarField::evaluate();
}


Foam::blackBodyRadiationFvPatchScalarField::
blackBodyRadiationFvPatchScalarField
(
    const blackBodyRadiationFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    emissivity_(ptf.emissivity_),
    kappa_(ptf.kappa_),
    Ta_(ptf.Ta_)
{}


Foam::blackBodyRadiationFvPatchScalarField::
blackBodyRadiationFvPatchScalarField
(
    const blackBodyRadiationFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    emissivity_(ptf.emissivity_),
    kappa_(ptf.kappa_),
    Ta_(ptf.Ta_)
{}


Foam::blackBodyRadiationFvPatchScalarField::
blackBodyRadiationFvPatchScalarField
(
    const blackBodyRadiationFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    emissivity_(ptf.emissivity_),
    kappa_(ptf.kappa_),
    Ta_(ptf.Ta_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::blackBodyRadiationFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //- Current patch temperature 
    scalarField Tp(*this);

    //- Internal field
    const scalarField Ti(patchInternalField());

    //- Inverse cell to patch face distances
    const scalarField invD = patch().deltaCoeffs();

    // Stefan-Boltzmann constant
    const scalar sigma_(constant::physicoChemical::sigma.value());

    //- Loop over all the faces of the boundary patch
    forAll(patch(), i)
    {   
        if (Ti[i] > 0)
        {
            //- Compute Fourier heat flux
            const scalar q = Foam::mag(-kappa_ * (Tp[i] - Ti[i]) * invD[i]);

            //- Compute outer cladding or metal/oxide interface temperature
            Tp[i] = Foam::pow(
                q/(sigma_ * emissivity_) + Foam::pow4(Ta_), 
                1.0/4.0
            );
        }
    }

    fvPatchField<scalar>::operator==(Tp);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::blackBodyRadiationFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("emissivity", emissivity_);
    os.writeEntry("kappa", kappa_);
    os.writeEntry("Ta", Ta_);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        blackBodyRadiationFvPatchScalarField
    );
}

// ************************************************************************* //
