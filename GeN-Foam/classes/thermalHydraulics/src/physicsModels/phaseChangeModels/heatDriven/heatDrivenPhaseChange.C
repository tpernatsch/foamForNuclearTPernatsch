/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "heatDrivenPhaseChange.H"
#include "FFPair.H"
#include "addToRunTimeSelectionTable.H"
#include "myOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeModels
{
    defineTypeNameAndDebug(heatDrivenPhaseChange, 0);
    addToRunTimeSelectionTable
    (
        phaseChangeModel,
        heatDrivenPhaseChange,
        phaseChangeModels
    );
}
}

const Foam::Enum
<
    Foam::phaseChangeModels::heatDrivenPhaseChange::mode
>
Foam::phaseChangeModels::heatDrivenPhaseChange::modeNames_
(
    {
        { 
            mode::conductionLimited, 
            "conductionLimited" 
        },
        { 
            mode::onePhaseDriven, 
            "onePhaseDriven" 
        },
        { 
            mode::twoPhaseDriven, 
            "twoPhaseDriven" 
        },
        { 
            mode::mixedDriven, 
            "mixedDriven" 
        }
    }
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModels::heatDrivenPhaseChange::heatDrivenPhaseChange
(
    FFPair& pair,
    const dictionary& dict
)
:
    phaseChangeModel
    (
        pair,
        dict
    ),
    mode_(modeNames_.get(this->get<word>("mode"))),
    drivingPhaseName_
    (
        (mode_ == heatDrivenPhaseChange::mode::onePhaseDriven) ?
        this->get<word>("drivingPhase") : ""
    )
{
    if (drivingPhaseName_ != "")
    {
        if  
        (
            drivingPhaseName_ != fluid1_.name()
        and drivingPhaseName_ != fluid2_.name()
        )
        {
            FatalErrorInFunction
            << "phaseChangeModel: phase " << drivingPhaseName_ << " not found!"
            << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseChangeModels::heatDrivenPhaseChange::correct() 
{
    //- Limit interfacial area so boiling can start
    //  (very crude, it's the best I have for now)
    limitInterfacialArea();

    //- Update saturation temperature
    correctInterfacialTemperature();

    //- Update latent heat
    myOps::storePrevIterIfRelax(L_);
    latentHeatPtr_->correctField(L_);
    L_.relax();
    
    //- Interfacial mass transfers for each side
    volScalarField dmdt1i
    (
        htc1_*iA_*(T1_-iT_)/L_
    );
    volScalarField dmdt2i
    (
        htc2_*iA_*(T2_-iT_)/L_
    );

    myOps::storePrevIterIfRelax(dmdt_);
    switch (mode_)
    {
        case heatDrivenPhaseChange::mode::conductionLimited : 
            dmdt_ = dmdt1i + dmdt2i;
            break;

        case heatDrivenPhaseChange::mode::twoPhaseDriven :
            dmdt_ = posPart(dmdt1i) + negPart(dmdt2i);
            break;

        case heatDrivenPhaseChange::mode::onePhaseDriven :
            dmdt_ = 
                (fluid1_.name() == drivingPhaseName_) ?
                dmdt1i : dmdt2i;      
            break;

        case heatDrivenPhaseChange::mode::mixedDriven :
            dmdt_ = 
                (fluid1_.name() == drivingPhaseName_) ?
                dmdt1i + negPart(dmdt2i) :
                posPart(dmdt1i) + dmdt2i;
            break;

        default : break;
    }
    limitMassTransfer();
    dmdt_.relax();
}


// ************************************************************************* //
