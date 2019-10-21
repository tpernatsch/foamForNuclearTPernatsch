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
    const FSPair& FSPair,
    const word region
)
:
    dragModel
    (
        objReg,
        dict,
        FSPair,
        region
    ),
    FSPair_(FSPair),
    Np_("numberOfPins", dimless, dict),
    Dp_("pinDiameter", dimLength, dict),
    Pp_("pinPitch", dimLength, dict),
    Dw_("wireDiameter", dimLength, dict),
    Lw_("wireLeadLen", dimLength, dict),
    wetWrapPer_("wetWrapPerimeter", dimLength, dict),
    wetPinPer_(Np_*constant::mathematical::pi*(Dp_+Dw_)),
    A_
    (
        wetPinPer_/(wetPinPer_+wetWrapPer_)
    )
{
    word principalAxisName(dict.get<word>("principalAxis"));
    if (principalAxisName == "localZ")
    {
        principalAxis_ = 2;
        transverseAxis1_ = 0;
        transverseAxis2_ = 1;
    }
    else if (principalAxisName == "localY")
    {
        principalAxis_ = 1;
        transverseAxis1_ = 0;
        transverseAxis2_ = 2;
    }
    else if (principalAxisName == "localX")
    {
        principalAxis_ = 0;
        transverseAxis1_ = 1;
        transverseAxis2_ = 2;
    }
    else
    {
        FatalErrorInFunction 
            << "Valid keywords for principalAxis are: localZ, localY, localX"
            << exit(FatalError);
    }

    dimensionedScalar B
    (
        sqrt(Pp_/Dp_) +
        pow
        (
            (
                7.6*(Dp_+Dw_)*sqr(Pp_/Dp_)
                /
                Lw_
            ),
            2.16
        )
    );
    B1_ = 64*sqrt(B);
    B2_ = 0.0816*pow(B, 0.9335);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volTensorField> Foam::dragModels::RehmeGunterShaw::Kd() const
{    
    tmp<volTensorField> tKd
    (
        new volTensorField
        (
            IOobject
            (
                "",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedTensor
            (
                "", 
                dimK, 
                tensor
                (
                    0, 0, 0,
                    0, 0, 0,
                    0, 0, 0
                )
            ),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volTensorField& Kd = tKd.ref();

    const volVectorField& lDh(FSPair_.structure().lDh());
    const volScalarField& magU(FSPair_.fluidRef().magU());
    tmp<volScalarField> pRe(mag(FSPair_.lRe().component(principalAxis_)));
    pRe = max(pRe, dimensionedScalar("", dimless, 10));
    tmp<volScalarField> tRe1(mag(FSPair_.lRe().component(transverseAxis1_)));
    tRe1 = max(tRe1, dimensionedScalar("", dimless, 10));
    tmp<volScalarField> tRe2(mag(FSPair_.lRe().component(transverseAxis2_)));
    tRe2 = max(tRe2, dimensionedScalar("", dimless, 10));

    volScalarField KdRehme
    (
        0.5*FSPair_.fluidRef().thermo().rho()*magU*
        (
            A_
            *
            (
                B1_/pRe()
                +
                B2_/pow(pRe(), 0.133)
            )
        )
        /
        lDh.component(principalAxis_)
    );
    Kd.replace(principalAxis_*4, KdRehme);

    volScalarField KdGunterShaw1
    (
        FSPair_.fluidRef().thermo().rho()*
        magU*
        0.96*pow(tRe1(),-0.145)
        /lDh.component(transverseAxis1_)
    );
    Kd.replace(transverseAxis1_*4, KdGunterShaw1);

    volScalarField KdGunterShaw2
    (
        FSPair_.fluidRef().thermo().rho()*
        magU*
        0.96*pow(tRe2(),-0.145)
        /lDh.component(transverseAxis2_)
    );
    Kd.replace(transverseAxis2_*4, KdGunterShaw2);
    
    Kd.correctBoundaryConditions();

    return tKd;
}


// ************************************************************************* //
