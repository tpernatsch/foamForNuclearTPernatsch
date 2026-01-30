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

#include "mappedDarcyVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "one.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedDarcyVelocityFvPatchVectorField::
mappedDarcyVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    isPorousRegion_()
{}


Foam::mappedDarcyVelocityFvPatchVectorField::
mappedDarcyVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    isPorousRegion_(dict.get<bool>("isPorousRegion"))
{

}


Foam::mappedDarcyVelocityFvPatchVectorField::
mappedDarcyVelocityFvPatchVectorField
(
    const mappedDarcyVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    isPorousRegion_(ptf.isPorousRegion_)
{}


Foam::mappedDarcyVelocityFvPatchVectorField::
mappedDarcyVelocityFvPatchVectorField
(
    const mappedDarcyVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    isPorousRegion_(ptf.isPorousRegion_)
{}


Foam::mappedDarcyVelocityFvPatchVectorField::
mappedDarcyVelocityFvPatchVectorField
(
    const mappedDarcyVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    isPorousRegion_(ptf.isPorousRegion_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedDarcyVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    
    // This patchfield needs to lookup the velocity from another region
    // and correct it so that the Darcy flux is conserved

    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
        this->patch().patch()
    );

    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch = refCast<const fvMesh>
    (
        nbrMesh
    ).boundary()[mpp.samplePolyPatch().index()];

    scalarField porosity =
        nbrPatch.lookupPatchField<volScalarField, scalar>("alpha");

    vectorField velocity =
        nbrPatch.lookupPatchField<volVectorField, vector>("U");

    

    vectorField darcyFlux(velocity);

    if(isPorousRegion_)
    {
        darcyFlux /= porosity;
    }
    else
    {
        darcyFlux *= porosity;
    }



    mpp.distribute(darcyFlux);

    this->operator==(darcyFlux);

    fixedValueFvPatchVectorField::updateCoeffs();



}


void Foam::mappedDarcyVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       mappedDarcyVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
