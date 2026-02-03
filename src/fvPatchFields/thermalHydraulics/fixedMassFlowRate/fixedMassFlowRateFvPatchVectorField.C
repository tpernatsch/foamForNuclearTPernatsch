/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2512                                                  |
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

#include "fixedMassFlowRateFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "one.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedMassFlowRateFvPatchVectorField::
fixedMassFlowRateFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    flowRate_(nullptr),
    rhoName_("thermo:rho"),
    rhoInlet_(0),
    volumetric_(false),
    extrapolateProfile_(false)
{}


Foam::fixedMassFlowRateFvPatchVectorField::
fixedMassFlowRateFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, IOobjectOption::NO_READ),
    flowRate_(nullptr),
    rhoName_("rho"),
    rhoInlet_(dict.getOrDefault<scalar>("rhoInlet", -VGREAT)),
    volumetric_(false),
    extrapolateProfile_
    (
        dict.getOrDefault<Switch>("extrapolateProfile", false)
    )
{
    flowRate_ =
        Function1<scalar>::NewIfPresent("volumetricFlowRate", dict, &db());

    if (flowRate_)
    {
        volumetric_ = true;
    }
    else
    {
        dict.readIfPresent("rho", rhoName_);
        flowRate_ =
            Function1<scalar>::NewIfPresent("massFlowRate", dict, &db());
    }

    if (!flowRate_)
    {
        FatalIOErrorInFunction(dict)
            << "Please supply either 'volumetricFlowRate' or"
            << " 'massFlowRate' (optional: with 'rho')" << nl
            << exit(FatalIOError);
    }

    // Value field required if mass based
    if (!this->readValueEntry(dict))
    {
        evaluate(Pstream::commsTypes::buffered);
    }
}


Foam::fixedMassFlowRateFvPatchVectorField::
fixedMassFlowRateFvPatchVectorField
(
    const fixedMassFlowRateFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    flowRate_(ptf.flowRate_.clone()),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    volumetric_(ptf.volumetric_),
    extrapolateProfile_(ptf.extrapolateProfile_)
{}


Foam::fixedMassFlowRateFvPatchVectorField::
fixedMassFlowRateFvPatchVectorField
(
    const fixedMassFlowRateFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    flowRate_(ptf.flowRate_.clone()),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    volumetric_(ptf.volumetric_),
    extrapolateProfile_(ptf.extrapolateProfile_)
{}


Foam::fixedMassFlowRateFvPatchVectorField::
fixedMassFlowRateFvPatchVectorField
(
    const fixedMassFlowRateFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    flowRate_(ptf.flowRate_.clone()),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_),
    volumetric_(ptf.volumetric_),
    extrapolateProfile_(ptf.extrapolateProfile_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class RhoType>
void Foam::fixedMassFlowRateFvPatchVectorField::updateValues
(
    const RhoType& rho
)
{
    const scalar t = db().time().timeOutputValue();

    const vectorField n(patch().nf());

    // Get alpha field

    const scalarField& alpha = 
        patch().lookupPatchField<volScalarField, scalar>("alpha");

    if (extrapolateProfile_)
    {
        vectorField Up(this->patchInternalField());

        // Patch normal extrapolated velocity
        scalarField nUp(n & Up);

        // Remove the normal component of the extrapolate patch velocity
        Up -= nUp*n;

        // Remove any reverse flow
        nUp = min(nUp, scalar(0));

        const scalar flowRate = flowRate_->value(t);
        const scalar estimatedFlowRate = -gSum(rho*(this->patch().magSf()*nUp));

        if (estimatedFlowRate > 0.5*flowRate)
        {
            nUp *= (mag(flowRate)/mag(estimatedFlowRate));
        }
        else
        {
            nUp -= ((flowRate - estimatedFlowRate)/gSum(rho*patch().magSf()));
        }

        // Add the corrected normal component of velocity to the patch velocity
        Up += nUp*n;
        Up /= alpha;

        // Correct the patch velocity
        this->operator==(Up);
    }
    else
    {
        const scalar avgU = -flowRate_->value(t)/gSum(alpha*rho*patch().magSf());
        operator==(avgU*n);
    }
}


void Foam::fixedMassFlowRateFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (volumetric_ || rhoName_ == "none")
    {
        updateValues(one{});
    }
    else
    {
        // Mass flow-rate
        if
        (
            const auto* rhop
          = patch().cfindPatchField<volScalarField>(rhoName_)
        )
        {
            updateValues(*rhop);
        }
        else
        {
            // Use constant density
            if (rhoInlet_ < 0)
            {
                FatalErrorInFunction
                    << "Did not find registered density field " << rhoName_
                    << " and no constant density 'rhoInlet' specified"
                    << exit(FatalError);
            }

            updateValues(rhoInlet_);
        }
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::fixedMassFlowRateFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    flowRate_->writeData(os);
    if (!volumetric_)
    {
        os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
        os.writeEntryIfDifferent<scalar>("rhoInlet", -VGREAT, rhoInlet_);
    }
    if (extrapolateProfile_)
    {
        os.writeEntry("extrapolateProfile", extrapolateProfile_);
    }
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       fixedMassFlowRateFvPatchVectorField
   );
}


// ************************************************************************* //
