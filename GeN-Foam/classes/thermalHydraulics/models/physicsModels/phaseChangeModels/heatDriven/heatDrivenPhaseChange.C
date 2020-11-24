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
#include "fluid.H"
#include "addToRunTimeSelectionTable.H"

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
    const dictionary& dict,
    const pimpleControl& pimple, 
    const fluid& fluid1,
    const fluid& fluid2,
    const volScalarField& p,
    const volScalarFieldPtrTable& htcs,
    volScalarField& dmdt,
    volScalarField& iT,
    volScalarField& iA
)
:
    phaseChangeModel
    (
        dict,
        pimple,
        fluid1,
        fluid2,
        p,
        htcs,
        dmdt,
        iT,
        iA
    ),
    mode_(modeNames_.get(this->get<word>("mode"))),
    drivingPhaseName_
    (
        (mode_ == heatDrivenPhaseChange::mode::onePhaseDriven) ?
        this->get<word>("drivingPhase") : ""
    ),
    L_
    (
        latentHeat_->L()
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

void Foam::phaseChangeModels::heatDrivenPhaseChange::correctMassTransfer() 
{
    //- Refs
    const volScalarField& T1(fluid1_.thermo().T());
    const volScalarField& T2(fluid2_.thermo().T());
    const volScalarField& H1i(*htcs_[fluid1_.name()+"."+fluid2_.name()]);
    const volScalarField& H2i(*htcs_[fluid2_.name()+"."+fluid1_.name()]);

    //- Limit interfacial area so boiling can start 
    //  (very crude, it's the best I have for now)
    this->limitInterfacialArea();

    //- Update latent heat
    latentHeat_->correct();
    
    volScalarField dmdt1i
    (
        H1i*iA_*(T1-iT_)/L_
    );
    volScalarField dmdt2i
    (
        H2i*iA_*(T2-iT_)/L_
    );

    dmdt_.storePrevIter();

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

    this->limitMassTransfer();

    dmdt_.relax();
}


// ************************************************************************* //
