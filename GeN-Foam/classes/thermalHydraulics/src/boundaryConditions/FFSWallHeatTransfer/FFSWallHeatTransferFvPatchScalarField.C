/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "FFSWallHeatTransferFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "rhoThermo.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FFSWallHeatTransferFvPatchScalarField::FFSWallHeatTransferFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    Tinf_(p.size(), 0.0),
    hWall_(p.size(), 0.0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::FFSWallHeatTransferFvPatchScalarField::FFSWallHeatTransferFvPatchScalarField
(
    const FFSWallHeatTransferFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    Tinf_(ptf.Tinf_, mapper),
    hWall_(ptf.hWall_, mapper)
{}


Foam::FFSWallHeatTransferFvPatchScalarField::FFSWallHeatTransferFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    Tinf_("Tinf", dict, p.size()),
    hWall_("hWall", dict, p.size())
{
    refValue() = Tinf_;
    refGrad() = 0.0;
    valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        evaluate();
    }
}


Foam::FFSWallHeatTransferFvPatchScalarField::FFSWallHeatTransferFvPatchScalarField
(
    const FFSWallHeatTransferFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    Tinf_(tppsf.Tinf_),
    hWall_(tppsf.hWall_)
{}


Foam::FFSWallHeatTransferFvPatchScalarField::FFSWallHeatTransferFvPatchScalarField
(
    const FFSWallHeatTransferFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    Tinf_(tppsf.Tinf_),
    hWall_(tppsf.hWall_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FFSWallHeatTransferFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
    Tinf_.autoMap(m);
    hWall_.autoMap(m);
}


void Foam::FFSWallHeatTransferFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const FFSWallHeatTransferFvPatchScalarField& tiptf =
        refCast<const FFSWallHeatTransferFvPatchScalarField>(ptf);

    Tinf_.rmap(tiptf.Tinf_, addr);
    hWall_.rmap(tiptf.hWall_, addr);
}


void Foam::FFSWallHeatTransferFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const phaseCompressibleTurbulenceModel& turbModel =
        db().lookupObject<phaseCompressibleTurbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                internalField().group()
            )
        );

    /*
    
    At the wall, hFluid = kFluid/cellHalfWidth == kFluid*deltaCoeffs()
    Instantaneous heat flux balance gives:
    Tpatch = 
            (1.0/(1.0+hFluid/hWall))*Tinf_ 
        +   (1 - (1.0/(1.0+hFluid/hWall)))*TpatchInternalField

    So:

        refValue = Tinf_
        valueFraction = (1.0/(1.0+hFluid/hWall))
        refGrad = 0;

    */

    valueFraction() =
        1.0/
        (
            1.0
          + turbModel.kappaEff(patch().index())*patch().deltaCoeffs()/hWall_
        );

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::FFSWallHeatTransferFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    Tinf_.writeEntry("Tinf", os);
    hWall_.writeEntry("hWall", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        FFSWallHeatTransferFvPatchScalarField
    );
}

// ************************************************************************* //
