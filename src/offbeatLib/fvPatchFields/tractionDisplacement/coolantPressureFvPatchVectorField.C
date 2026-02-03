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

#include "coolantPressureFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

#include "globalOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coolantPressureFvPatchVectorField::
coolantPressureFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(p, iF),
    coolantPressure_(p.size(), 0.0)
{}


coolantPressureFvPatchVectorField::
coolantPressureFvPatchVectorField
(
    const coolantPressureFvPatchVectorField& tvtdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    tractionDisplacementFvPatchVectorField(tvtdpvf, p, iF, mapper),
    coolantPressure_(tvtdpvf.coolantPressure_)
{}


coolantPressureFvPatchVectorField::
coolantPressureFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    tractionDisplacementFvPatchVectorField(p, iF, dict, false),
    coolantPressure_(p.size(), 0)
{
    // Read fluid pressure (either list or fixed pressure)
    if(dict.found("coolantPressureList"))
    {

#ifdef OPENFOAMFOUNDATION            
            coolantPressureList_.set
            ( 
                new Function1s::Table<scalar>
                (
                    "coolantPressureList", dict.subDict("coolantPressureList")
                )
            );
#elif OPENFOAMESI
            word listCoeffsName = "coolantPressureList";
            if(dict.found("coolantPressureListCoeffs"))
            {
                listCoeffsName = "coolantPressureListCoeffs";
            }
            coolantPressureList_.reset
            ( 
                new scalarTable
                (
                    "coolantPressureList", dict.subDict(listCoeffsName)
                )
            );
#endif

        // Set pressure value
        coolantPressure_ = 
        coolantPressureList_->value
        (
            patch().boundaryMesh().mesh().time().timeToUserTime
            (
                this->db().time().value()
            )
        );
    }
    else if (dict.found("coolantPressure"))
    {
        coolantPressure_ = scalarField("coolantPressure", dict, p.size());
    }
    else
    {
        FatalErrorInFunction() 
        << "coolantPressureFvPatchVectorField BC for " << patch().name() << " requires:" << nl
        << "- \"coolantPressureList\" or" << nl
        << "- \"coolantPressure\"" << abort(FatalError) << endl;
    }

    // Read initial pressure value if present (for accurate restart)
    if(dict.found("pressure"))
    {
        pressure_ = scalarField("pressure", dict, p.size());
    }
    else
    {
        pressure_ = coolantPressure_;
    }
}


coolantPressureFvPatchVectorField::
coolantPressureFvPatchVectorField
(
    const coolantPressureFvPatchVectorField& tvtdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(tvtdpvf, iF),
    coolantPressure_(tvtdpvf.coolantPressure_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void coolantPressureFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    tractionDisplacementFvPatchVectorField::autoMap(m);
}


void coolantPressureFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    tractionDisplacementFvPatchVectorField::rmap(ptf, addr);
}


void coolantPressureFvPatchVectorField::updateTraction()
{
    // Fluid pressure
    if(coolantPressureList_.valid())
    {
        coolantPressure_ = 
        coolantPressureList_->value
        (
            this->db().time().timeToUserTime
            (
                this->db().time().value()
            )
        );
    }

    // Update total pressure
    pressure_ = coolantPressure_;
}


void coolantPressureFvPatchVectorField::write(Ostream& os) const
{
    tractionDisplacementFvPatchVectorField::write(os);

    if (coolantPressureList_.valid())
    {
        // os.writeKeyword("coolantPressureList") << nl;
        // os << token::BEGIN_BLOCK << nl;
#ifdef OPENFOAMFOUNDATION  
        coolantPressureList_->write(os);
#elif OPENFOAMESI
        coolantPressureList_->writeData(os);
#endif
        // os << token::END_BLOCK << nl;
    }
    else
    {
#ifdef OPENFOAMFOUNDATION  
        writeEntry<scalarField>(os, "coolantPressure", coolantPressure_);
#elif OPENFOAMESI
        coolantPressure_.writeEntry( "coolantPressure", os);
#endif 
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    coolantPressureFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
