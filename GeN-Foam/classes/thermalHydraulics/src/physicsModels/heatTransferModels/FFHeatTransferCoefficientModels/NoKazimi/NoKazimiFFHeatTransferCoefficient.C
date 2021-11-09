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

#include "FFPair.H"
#include "NoKazimiFFHeatTransferCoefficient.H"
#include "mathematicalConstants.H"
#include "physicoChemicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FFHeatTransferCoefficientModels
{
    defineTypeNameAndDebug(NoKazimi, 0);
    addToRunTimeSelectionTable
    (
        FFHeatTransferCoefficientModel, 
        NoKazimi, 
        FFHeatTransferCoefficientModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FFHeatTransferCoefficientModels::NoKazimi::NoKazimi
(
    const FFPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    FFHeatTransferCoefficientModel
    (
        pair,
        dict,
        objReg
    ),
    alpha_(bulkFluid_.normalized()),
    rho_(bulkFluid_.rho()),
    kappa_(bulkFluid_.kappa()),
    Cp_(bulkFluid_.Cp()),
    Pr_(bulkFluid_.Cp()),
    magU_(bulkFluid_.magU()),
    p_(pair.mesh().lookupObject<volScalarField>("p")),
    iA_(pair.iA()),
    iT_(pair.iT()),
    coeff_
    (
        sqrt
        (
            // W/1000 g/mol->kg/mol
            bulkFluid_.thermo().W()().average().value()/1000.0/ 
            (
                2.0*constant::mathematical::pi*
                constant::physicoChemical::R.value()
            )
        )
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::FFHeatTransferCoefficientModels::NoKazimi::value
(
    const label& celli
) const
{
    scalar alpha(max(alpha_[celli], 1e-9));
    const scalar& rho(rho_[celli]);
    const scalar& L(pair_.L()[celli]);
    return 
        min
        (
            iA_[celli]*kappa_[celli]/(alpha*(1.0-alpha))
        +   27.73*bulkFluid_[celli]*rho*magU_[celli]*Cp_[celli]*Pr_[celli],
            coeff_*sqr(rho*L)/(p_[celli]*sqrt(iT_[celli]))
        );
}


// ************************************************************************* //
