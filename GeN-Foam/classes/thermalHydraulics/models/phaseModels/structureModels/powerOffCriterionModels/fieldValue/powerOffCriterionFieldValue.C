/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "powerOffCriterionFieldValue.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace powerOffCriterionModels
{
    defineTypeNameAndDebug(fieldValue, 0);
    addToRunTimeSelectionTable
    (
        powerOffCriterionModel, 
        fieldValue, 
        powerOffCriterionModels
    );
}
}

const Foam::Enum
<
    Foam::powerOffCriterionModels::fieldValue::fieldOp
>
Foam::powerOffCriterionModels::fieldValue::fieldOpNames_
(
    {
        { 
            fieldOp::max,
            "max" 
        },
        { 
            fieldOp::min,
            "min" 
        }
    }
);

const Foam::Enum
<
    Foam::powerOffCriterionModels::fieldValue::criterion
>
Foam::powerOffCriterionModels::fieldValue::criterionNames_
(
    {
        { 
            criterion::above, 
            "valueAboveThreshold" 
        },
        { 
            criterion::below, 
            "valueBelowThreshold" 
        }
    }
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerOffCriterionModels::fieldValue::fieldValue
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    powerOffCriterionModel
    (
        mesh,
        dict
    ),
    fieldName_(this->get<word>("fieldName")),
    fieldOp_
    (
        fieldOpNames_.get
        (
            this->get<word>
            (
                "fieldOperation"
            )
        )
    ),
    criterion_
    (
        criterionNames_.get
        (
            this->get<word>
            (
                "criterion"
            )
        )
    ),
    threshold_(this->get<scalar>("threshold")),
    timeDelay_(this->lookupOrDefault<scalar>("timeDelay", 0.0)),
    time0_(0.0),
    fieldPtr_(nullptr),
    InfoFlag_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerOffCriterionModels::fieldValue::~fieldValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::powerOffCriterionModels::fieldValue::powerOffCriterion()
{   
    if (fieldPtr_ == nullptr)
    {
        fieldPtr_ = &mesh_.lookupObject<volScalarField>(fieldName_);
    }

    scalar value(0.0);
    bool flag(false);

    switch(fieldOp_)
    {
        case fieldOp::max :
        {
            value = Foam::max(*fieldPtr_).value();
            break;
        }
        case fieldOp::min :
        {
            value = Foam::min(*fieldPtr_).value();
            break;
        }
    }

    switch(criterion_)
    {
        case criterion::above :
        {
            flag = (value >= threshold_);
            break;
        }
        case criterion::below :
        {
            flag = (value <= threshold_);
            break;
        }
    }

    if (flag)
    {
        if (timeDelay_ != 0.0)
        {
            scalar time(mesh_.time().timeOutputValue());
            if (time0_ == 0.0) time0_ = time;
            flag = (time >= (time0_ + timeDelay_));
        }
        else if (!InfoFlag_) time0_ = mesh_.time().timeOutputValue();
    }
    if (flag and !InfoFlag_) InfoFlag_ = true;
    if (InfoFlag_) 
            Info << "Power off at t = " << time0_+timeDelay_ << " s" << endl;
    return flag;
}


// ************************************************************************* //
