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

#include "plenumSpringPressureFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

#include "gapGasModel.H"
#include "globalOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plenumSpringPressureFvPatchVectorField::
plenumSpringPressureFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(p, iF),
    plenumSpringBase(p)
{}


plenumSpringPressureFvPatchVectorField::
plenumSpringPressureFvPatchVectorField
(
    const plenumSpringPressureFvPatchVectorField& tvtdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    tractionDisplacementFvPatchVectorField(tvtdpvf, p, iF, mapper),
    plenumSpringBase(p, *this)
{}


plenumSpringPressureFvPatchVectorField::
plenumSpringPressureFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    tractionDisplacementFvPatchVectorField(p, iF, dict, false),
    plenumSpringBase(p, dict)
{
    // Read initial pressure value if present (for accurate restart)
    if(dict.found("pressure"))
    {
        pressure_ = scalarField("pressure", dict, p.size());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void plenumSpringPressureFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    tractionDisplacementFvPatchVectorField::autoMap(m);
}


void plenumSpringPressureFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    tractionDisplacementFvPatchVectorField::rmap(ptf, addr);
}


void plenumSpringPressureFvPatchVectorField::updateTraction()
{
    const globalOptions& globalOpt
    (db().lookupObject<globalOptions>("globalOptions"));

    scalar angularFraction = globalOpt.angularFraction();

    // Include the gap gas pressure
    const gapGasModel& gapGas
    = db().lookupObject<gapGasModel>("gapGas");    

    // Calculate total spring elongation (positive if spring contracts)
    scalar Dtot(plenumSpringBase::springElongation());

    scalar Atot = gSum(patch().magSf());

    scalarField springPressure
    (patch().size(), plenumSpringBase::springModulus()*Dtot/(Atot/angularFraction));

    // Update total pressure
    pressure_ = springPressure + gapGas.p();  
}


void plenumSpringPressureFvPatchVectorField::write(Ostream& os) const
{
    tractionDisplacementFvPatchVectorField::write(os);
    plenumSpringBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    plenumSpringPressureFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
