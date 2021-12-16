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
#include "LookUpTableCHF.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace CriticalHeatFluxModels
{
    defineTypeNameAndDebug(LookUpTableCHF, 0);
    addToRunTimeSelectionTable
    (
        CHFModel,
        LookUpTableCHF,
        CriticalHeatFluxModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CriticalHeatFluxModels::LookUpTableCHF::LookUpTableCHF
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    CHFModel
    (
        pair,
        dict,
        objReg
    ),
    quality_(pair.fluidRef().flowQuality()),
    p_(pair.mesh().lookupObject<volScalarField>("p")),
    rhoL_(pair.mesh().lookupObject<volScalarField>("thermo:rho.liquid")),
    rhoV_(pair.mesh().lookupObject<volScalarField>("thermo:rho.vapour")),
    uL_(pair.fluidRef().magU()),
    uV_(pair.mesh().lookupObject<volScalarField>("magU.vapour")),
    normalizedL_(pair.fluidRef().normalized())
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::CriticalHeatFluxModels::LookUpTableCHF::value
(
    const label& celli
) const
{
    const scalar& xi(quality_[celli]);
    const scalar& pi(p_[celli]);
    const scalar& rhoLi(rhoL_[celli]);
    const scalar& rhoVi(rhoV_[celli]);
    const scalar& uLi(uL_[celli]);
    const scalar& uVi(uV_[celli]);
    const scalar& aLi(normalizedL_[celli]);


    scalar massFlowi(rhoLi*uLi*aLi+rhoVi*uVi*(1-aLi));

    // Temporary Gauthier 




    return 0;
}

// ************************************************************************* //
