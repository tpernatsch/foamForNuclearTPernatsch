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

#include "FSPair.H"
#include "CachardVapourFSHeatTransferCoefficient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace FSHeatTransferCoefficientModels
{
    defineTypeNameAndDebug(CachardVapour, 0);
    addToRunTimeSelectionTable
    (
        FSHeatTransferCoefficientModel, 
        CachardVapour, 
        FSHeatTransferCoefficientModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSHeatTransferCoefficientModels::CachardVapour::CachardVapour
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    FSHeatTransferCoefficientModel
    (
        pair,
        dict,
        objReg
    ),
    p_(pair.mesh().lookupObject<volScalarField>("p")),
    alpha_(pair.fluidRef().normalized()),
    kg_(pair.mesh().lookupObject<volScalarField>("kappa.vapour")),
    Dh_(pair.fluidRef().Dh()),
    gravity_(9.81)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
IMPLEMENTATION NOTES


*/

Foam::scalar Foam::FSHeatTransferCoefficientModels::CachardVapour::value
(
    const label& celli
) const
{
    const scalar& kgi(kg_[celli]);
    const scalar& Dhi(Dh_[celli]);
    const scalar DRi = 0.01;
    const scalar& alphai(alpha_[celli]);
    const scalar& pi(p_[celli]);


    // Film Thickness delta 
    scalar deltai(Dhi/2*(pow(1+alphai*((4/Foam::constant::mathematical::pi)*pow(pi/max(DRi,1e-6),2)-1),0.5)-1));
    scalar hCachard = 2*kgi/max(deltai,1e-6);

    return hCachard;
}

// ************************************************************************* //
