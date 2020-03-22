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

#include "ChurchillGunterShaw.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(ChurchillGunterShaw, 0);
    addToRunTimeSelectionTable(dragModel, ChurchillGunterShaw, FSDragModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::ChurchillGunterShaw::ChurchillGunterShaw
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
    surfaceRoughness_(this->get<scalar>("surfaceRoughness"))
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
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dragModels::ChurchillGunterShaw::correctKd(volTensorField& Kd) const
{    
    //- Init refs
    const volVectorField& lRe(FSPair_->lRe());
    const volVectorField& lDh(FSPair_->structure().lDh());
    const volScalarField& magU(FSPair_->fluidRef().magU());
    volScalarField halfRhoMagU(0.5*FSPair_->fluidRef().thermo().rho()*magU);

    //- Direction labeles
    label d0(lcmpt_[0]); //- Transverse direction 0
    label d1(lcmpt_[1]); //- Transverse direction 1
    label d2(lcmpt_[2]); //- Principal direction

    //- Set field cell-by-cell. It is appreciably faster than constructing
    //  whole-mesh fields to perform region-by-region operations
    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        const vector& lRei(lRe[celli]);
        const vector& lDhi(lDh[celli]);
        const scalar& halfRhoMagUi(halfRhoMagU[celli]);
        const scalar& Re0(lRei[d0]);
        const scalar& Re1(lRei[d1]);
        const scalar& Re2(lRei[d2]);
        tensor& Kdi(Kd[celli]);

        //- Churchill pressure drop for principal direction (i.e. d2). What
        //  about the weird indexing? If the principal direction is localX,
        //  then d2 = 0 and thus I am setting the xx component of the local
        //  drag tensor. If the principal direction is localY, then d2 = 1
        //  and thus I am setting the yy component of the local drag tensor
        //  (i.e. component 4). If the principal direction is localZ, then
        //  d2 = 2 and thus I am setting the zz component of the local drag
        //  tensor (i.e. component 8). Here is a reminder on tensor component
        //  numbering
        /*
            | xx xy xz |   | 0 1 2 |
            | yx yy yz | = | 3 4 5 |
            | zx zy zz |   | 6 7 8 |
        */

        //- Churchill coeffs
        scalar A
        (
            pow
            (
                -2.457*
                Foam::log
                (
                    pow(7.0/Re2, 0.9)
                +   0.27*surfaceRoughness_
                ), 
                16
            ) 
        );
        scalar B(pow(37530.0/Re2, 16));

        Kdi[d2*4] =
        (
            8.0*halfRhoMagUi*
            pow
            (
                pow(8.0/Re2, 12)
            +   1.0/pow(A+B, 1.5),
                1.0/12
            )/
            lDhi[d2]
        );

        //- Gunter-Shaw pressure drop for non-principal direction, d0 and d1
        Kdi[d0*4] = halfRhoMagUi*0.96*pow(Re0,-0.145)/lDhi[d0];
        Kdi[d1*4] = halfRhoMagUi*0.96*pow(Re1,-0.145)/lDhi[d1];
    }
}


// ************************************************************************* //
