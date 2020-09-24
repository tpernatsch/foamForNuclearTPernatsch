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

#include "NoKazimiFFHeatTransfer.H"
#include "mathematicalConstants.H"
#include "physicoChemicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(NoKazimiFFHeatTransfer, 0);
    addToRunTimeSelectionTable
    (
        heatTransferModel, 
        NoKazimiFFHeatTransfer, 
        FFHeatTransferModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::NoKazimiFFHeatTransfer::
NoKazimiFFHeatTransfer
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FFPair& FFPair,
    const word nameBulk,
    const word nameInterface
)
:
    heatTransferModel
    (
        objReg,
        dict,
        FFPair,
        nameBulk,
        nameInterface
    ),
    bulkFluid_
    (
        (
            FFPair.fluid1().name() == nameBulk
        ) ?
        FFPair.fluid1() :
        FFPair.fluid2()
    ),
    otherFluid_
    (
        (
            FFPair.fluid1().name() == nameBulk
        ) ?
        FFPair.fluid2() :
        FFPair.fluid1()
    ),
    p_(FFPair.fluid1().mesh().lookupObject<volScalarField>("p")),
    iA_
    (
        FFPair.fluid1().mesh().lookupObject<volScalarField>
        (
            "iA."+FFPair.fluid1().name()+"."+FFPair.fluid2().name()
        )
    ),
    iT_(FFPair.fluid1().mesh().lookupObject<volScalarField>("T.interface")),
    coeff_
    (
        sqrt
        (
            bulkFluid_.thermo().W()().average()/1000.0/ // W/1000 g/mol->kg/mol
            (
                2.0*constant::mathematical::pi*constant::physicoChemical::R
            )
        )*dimensionedScalar
        (
            "",
            dimArea/dimTime/dimTime/dimTemperature,
            1.0
        )
    )
{
    if 
    (
            FFPair.fluid1().name() != nameBulk 
        and FFPair.fluid2().name() != nameBulk
    )
    {
        FatalErrorInFunction
            << "Phase " << nameBulk << " not found!"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::heatTransferModels::NoKazimiFFHeatTransfer::correctHtc
(
    volScalarField& htc
) const
{    
    //- Refs
    const volScalarField& alpha(bulkFluid_.normalized());
    const volScalarField& rho(bulkFluid_.thermo().rho());
    const volVectorField& U(bulkFluid_.U());
    
    //- Tmps to which I cannot take refs but need to be evaulated on the fly
    tmp<volScalarField> tkappa(bulkFluid_.thermo().kappa());
    volScalarField& kappa = tkappa.ref();
    tmp<volScalarField> tCp(bulkFluid_.thermo().Cp());
    volScalarField& Cp = tCp.ref();
    
    //- Latent heat
    volScalarField H(mag(bulkFluid_.thermo().hc()-otherFluid_.thermo().hc()));

    volScalarField htcA
    (
        iA_*kappa/max(alpha*(1.0-alpha), 1e-2)
    +   27.73*bulkFluid_*rho*mag(U)*Cp*bulkFluid_.Pr()
    );

    volScalarField htcB
    (
        coeff_*sqr(rho*H)/(p_*sqrt(iT_))
    );

    htc = min(htcA, htcB);
}


// ************************************************************************* //
