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
    correctLatentHeat_
    (
        this->lookupOrDefault<bool>("correctLatentHeat", false)
    ),
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

void Foam::phaseChangeModels::heatDrivenPhaseChange::correctMassTransfer() 
{
    //- Refs
    const volScalarField& he1(fluid1_.thermo().he());
    const volScalarField& he2(fluid2_.thermo().he());
    const volScalarField& T1(fluid1_.thermo().T());
    const volScalarField& T2(fluid2_.thermo().T());
    const volScalarField& H1i(*htcs_[fluid1_.name()+"."+fluid2_.name()]);
    const volScalarField& H2i(*htcs_[fluid2_.name()+"."+fluid1_.name()]);

    //- Limit interfacial area so boiling can start 
    //  (very crude, it's the best I have for now)
    this->limitInterfacialArea();

    //- correctLatentHeat_ =
    //  compute the latent heat of vaporization in a manner consistent
    //  with how massTransfer related enthalpy is added to the energy
    //  equations. In particular, assuming 2 is vapour and 1 is liquid
    //  (but it works the other way around too):
    //  * when evaporating, mass is removed from the liquid at its current
    //    enthalpy and is added to the vapour at vapour saturation
    //    entalphy;
    //  * when condensing, mass is removed from the vapour at its current
    //    enthalpy and is added to the liquid at liquid saturation
    //    enthalpy.
    //    For example, in the case of condensing super-heated steam, its
    //    enthalpy will be larger than saturation enthalpy. However, if
    //    you remove mass from the steam at an enthalpy lower than its
    //    current enthalpy, the remaining steam will inevitably be
    //    hotter.
    //  The un-adjusted latent heat itself would be
    //  fluid2_.thermo().hc() - fluid1_.thermo().hc(), with hc() being
    //  the enthalpy of formation specified under the Hf keyword in the
    //  thermophysical dict of each phase

    //- Please note that the code here does not make use of the isLiquid or
    //  isGas methods of the fluids to determine which is the liquid and 
    //  which is the vapour. It does that via the sign of the latent heat
    volScalarField L(fluid2_.thermo().hc() - fluid1_.thermo().hc());
    if (correctLatentHeat_)
    {
        L += 
            neg0(dmdt_)*(he2 - fluid2_.thermo().he(p_, iT_))
        -   pos0(dmdt_)*(he1 - fluid1_.thermo().he(p_, iT_));

        Info<< "L (avg min max) ="
        << " " << L.weightedAverage(mesh_.V()).value()
        << " " << min(L).value()
        << " " << max(L).value()
        << " J/kg" << endl;
    }

    scalar f(this->relaxationFactor());
    
    volScalarField dmdt1i
    (
        pos0(fluid1_-fluid1_.residualAlpha())*
        H1i*iA_*(T1-iT_)/L
    );
    volScalarField dmdt2i
    (
        pos0(fluid2_-fluid2_.residualAlpha())*
        H2i*iA_*(T2-iT_)/L
    );

    switch (mode_)
    {
        case heatDrivenPhaseChange::mode::conductionLimited : 

            dmdt_ = 
                (1.0-f)*dmdt_ + f*(dmdt1i + dmdt2i);
            break;

        case heatDrivenPhaseChange::mode::twoPhaseDriven :
            dmdt_ = 
                (1.0-f)*dmdt_ + f*(posPart(dmdt1i) + negPart(dmdt2i));
            break;

        case heatDrivenPhaseChange::mode::onePhaseDriven :
            dmdt_ = 
                (1.0-f)*dmdt_ 
            +   
                f*
                (
                    (fluid1_.name() == drivingPhaseName_) ?
                    dmdt1i 
                    :
                    dmdt2i 
                );      
            break;

        default : break;
    }

    this->limitMassTransfer();
}


// ************************************************************************* //
