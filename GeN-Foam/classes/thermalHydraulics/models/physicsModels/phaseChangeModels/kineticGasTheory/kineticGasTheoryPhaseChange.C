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

#include "kineticGasTheoryPhaseChange.H"
#include "fluid.H"
#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "physicoChemicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeModels
{
    defineTypeNameAndDebug(kineticGasTheoryPhaseChange, 0);
    addToRunTimeSelectionTable
    (
        phaseChangeModel,
        kineticGasTheoryPhaseChange,
        phaseChangeModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModels::kineticGasTheoryPhaseChange::kineticGasTheoryPhaseChange
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
    sigma_("sigma", dimless, 1.0), //this),
    lambdaEvap_("evaporationCoeff", dimless, this),
    lambdaCond_("condensationCoeff", dimless, this),
    coeff_
    (
        2.0*sigma_/(2.0-sigma_)*
        sqrt
        (
            fluid1_.thermo().W()().average()/1000.0/ // W/1000 g/mol->kg/mol
            (
                2.0*constant::mathematical::pi*constant::physicoChemical::R
            )
        )
    ),
    //alphaDryout_("alphaDryout", dimless, this),
    Tl_
    (
        (fluid1_.isLiquid()) ? fluid1_.thermo().T() : fluid2_.thermo().T()
    ),
    Tv_
    (
        (fluid1_.isGas()) ? fluid1_.thermo().T() : fluid2_.thermo().T()
    ),
    alphal_
    (
        (fluid1_.isLiquid()) ? fluid1_ : fluid2_
    ),
    alphav_
    (
        (fluid1_.isGas()) ? fluid1_ : fluid2_
    )
{
    dimensionedScalar M1(fluid1_.thermo().W()().average().value());
    dimensionedScalar M2(fluid2_.thermo().W()().average().value());

    if (M1.value() != M2.value())
    {
        FatalErrorInFunction
        << "phaseChangeModel: inconsistent molar masses between "
        << fluid1_.name() << " and " << fluid2_.name() << exit(FatalError); 
    }

    if (sigma_.value() > 1.0 or sigma_.value() < 0.0)
    {
        FatalErrorInFunction
        << "phaseChangeModel: sigma must be between 0 and 1" 
        << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseChangeModels::kineticGasTheoryPhaseChange::correctMassTransfer() 
{
    //- Limit interfacial area so boiling can start 
    //  (very crude, it's the best I have for now)
    this->limitInterfacialArea();

    scalar f(this->relaxationFactor());

    //volScalarField alphaDiff(alphav/(alphav+alphal)-alphaDryout_);

    //- This should be implict w.r.t. alphav, Tl, Tv, but yeah, good luck
    //  trying to do that... Maybe I can get around the alpha limitation
    //  (e.g. by copying alphav and solving it with a dmdt_ expressed in terms
    //  of an Sp(G, alphav) to get at estimate for alphav at the next 
    //  iteration), but definitely not a priority now. Thus, for the time
    //  being, the alphas that appear in the equation are limited so that
    //  G can be different than 0 after conditions for phase change are met
    volScalarField G
    (
        coeff_*max(alphal_, 1e-2)*max(alphav_, 1e-2)*
        (
            lambdaEvap_*posPart(saturation_->pSat(Tl_)-p_) 
        +   lambdaCond_*negPart(saturation_->pSat(Tv_)-p_)
        )/sqrt(iT_)
    );

    dmdt_ = 
        (1.0-f)*dmdt_ + f*iA_*G;    //- fiGA_ would have been nicer for 
                                    //  obvious italian reasons

    this->limitMassTransfer();
}


// ************************************************************************* //
