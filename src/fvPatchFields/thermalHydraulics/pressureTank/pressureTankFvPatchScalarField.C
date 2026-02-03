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

#include "pressureTankFvPatchScalarField.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureTankFvPatchScalarField::pressureTankFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    p0_(0),
    level_(0),
    area_(1),
    currTimeIndex_(0),
    totalMass_(0),
    deltaM_(0)
{}
 
 
Foam::pressureTankFvPatchScalarField::pressureTankFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    p0_(dict.get<scalar>("p0")),
    level_(dict.get<scalar>("level")),
    area_(dict.get<scalar>("area")),
    currTimeIndex_(0),
    totalMass_(0),
    deltaM_(0)
{}

 
 
Foam::pressureTankFvPatchScalarField::pressureTankFvPatchScalarField
(
    const pressureTankFvPatchScalarField& sppsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(sppsf, p, iF, mapper),
    p0_(sppsf.p0_),
    level_(sppsf.level_),
    area_(sppsf.area_),
    currTimeIndex_(sppsf.currTimeIndex_),
    totalMass_(sppsf.totalMass_),
    deltaM_(sppsf.deltaM_)
{}
 
 
Foam::pressureTankFvPatchScalarField::pressureTankFvPatchScalarField
(
    const pressureTankFvPatchScalarField& sppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(sppsf, iF),
    p0_(sppsf.p0_),
    level_(sppsf.level_),
    area_(sppsf.area_),
    currTimeIndex_(sppsf.currTimeIndex_),
    totalMass_(sppsf.totalMass_),
    deltaM_(sppsf.deltaM_)
{}
 
 
Foam::pressureTankFvPatchScalarField::pressureTankFvPatchScalarField
(
    const pressureTankFvPatchScalarField& sppsf
)
:
    fixedValueFvPatchScalarField(sppsf),
    p0_(sppsf.p0_),
    level_(sppsf.level_),
    area_(sppsf.area_),
    currTimeIndex_(sppsf.currTimeIndex_),
    totalMass_(sppsf.totalMass_),
    deltaM_(sppsf.deltaM_)

{}
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pressureTankFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalarField& density =
        patch().lookupPatchField<volScalarField, scalar>("thermo:rho"); 

    if (db().time().timeIndex() == 1)
        totalMass_ = average(density)* level_ * area_;


    if (currTimeIndex_ != db().time().timeIndex())
    {
        totalMass_ = totalMass_+deltaM_;
        currTimeIndex_ = db().time().timeIndex();
        Info << "Total mass in tank " << patch().name()<< ": " << totalMass_ << endl;
    }

    scalar t = db().time().value();
    scalar deltaT = db().time().deltaTValue();


    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>("alphaRhoPhi");
    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);
    


    deltaM_ = deltaT*sum(phip);

    level_ = (totalMass_+deltaM_)/(average(density)*area_);
   
   
    uniformDimensionedVectorField g = this->db().lookupObject<uniformDimensionedVectorField>("g");

    const surfaceScalarField& ghf =
        db().lookupObject<surfaceScalarField>("ghf");

    const fvsPatchField<scalar>& ghfp =
        patch().patchField<surfaceScalarField, scalar>(ghf);


    scalar newPressure(0);
    
    if(level_>0)
        newPressure= p0_+average(density)*(mag(g.value()))*(level_-average(ghfp));
    else    
        newPressure = p0_;

    operator==(newPressure);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::pressureTankFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    os.writeEntry("level", level_);
    os.writeEntry("area", area_);    
    os.writeEntry("p0", p0_);
    os.writeEntry("totalMass", totalMass_);

    fvPatchField<scalar>::writeValueEntry(os);
}


namespace Foam 
{ 
    makePatchTypeField(fvPatchScalarField, pressureTankFvPatchScalarField); 
}

// ************************************************************************* //
