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

#include "RehmeGunterShaw.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(RehmeGunterShaw, 0);
    addToRunTimeSelectionTable(dragModel, RehmeGunterShaw, FSDragModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::RehmeGunterShaw::RehmeGunterShaw
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FSPair& FSPair
)
:
    dragModel
    (
        objReg,
        dict,
        FSPair
    ),
    lcmpt_(vector::zero),
    Np_(this->get<label>("numberOfPins")),
    Dp_(this->get<scalar>("pinDiameter")),
    Pp_(this->get<scalar>("pinPitch")),
    Dw_(this->get<scalar>("wireDiameter")),
    Lw_(this->get<scalar>("wireLeadLen")),
    wetWrapPer_(this->get<scalar>("wetWrapPerimeter")),
    wetPinPer_(Np_*constant::mathematical::pi*(Dp_+Dw_)),
    A_
    (
        wetPinPer_/(wetPinPer_+wetWrapPer_)
    )
{
    word principalAxis(this->get<word>("principalAxis"));
    if (principalAxis == "localX")
    {
        lcmpt_ = vector(1,2,0);
    }
    else if (principalAxis == "localY")
    {
        lcmpt_ = vector(2,0,1);
    }
    else if (principalAxis == "localZ")
    {
        lcmpt_ = vector(0,1,2);
    }
    else
    {
        FatalErrorInFunction 
            << "Valid keywords for principalAxis are: localZ, localY, localX"
            << exit(FatalError);
    }

    scalar B(sqrt(Pp_/Dp_) + pow((7.6*(Dp_+Dw_)*sqr(Pp_/Dp_)/Lw_), 2.16));
    B1_ = 64*sqrt(B);
    B2_ = 0.0816*pow(B, 0.9335);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dragModels::RehmeGunterShaw::correctKd(volTensorField& Kd) const
{    
    //- Init refs
    const volScalarField& alpha(FSPair_->fluidRef());
    const volVectorField& lRe(FSPair_->lRe());
    const volVectorField& lDh(FSPair_->structure().lDh());
    const volScalarField& magU(FSPair_->fluidRef().magU());
    const volScalarField& rho(FSPair_->fluidRef().thermo().rho());

    //- Direction labeles
    label d0(lcmpt_[0]); //- Transverse direction 0
    label d1(lcmpt_[1]); //- Transverse direction 1
    label d2(lcmpt_[2]); //- Principal direction

    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        const vector& lRei(lRe[celli]);
        const vector& lDhi(lDh[celli]);
        const scalar& Re0(lRei[d0]);
        const scalar& Re1(lRei[d1]);
        const scalar& Re2(lRei[d2]);
        tensor& Kdi(Kd[celli]);
        
        scalar alphaRhoMagUi(0.5*alpha[celli]*rho[celli]*magU[celli]);

        //- Rehme correlation in principalAxis local dir (i.e. d2)
        Kdi[d2*4] = alphaRhoMagUi*(A_*(B1_/Re2 + B2_/pow(Re2, 0.133)));

        //- Gunter-Shaw for transverse axes
        Kdi[d0*4] = alphaRhoMagUi*0.96*pow(Re0, -0.145)/lDhi[d0];
        Kdi[d1*4] = alphaRhoMagUi*0.96*pow(Re1, -0.145)/lDhi[d1];
    }
}


// ************************************************************************* //
