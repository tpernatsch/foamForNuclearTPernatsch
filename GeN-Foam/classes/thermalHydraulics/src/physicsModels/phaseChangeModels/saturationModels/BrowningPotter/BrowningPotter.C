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

#include "FFPair.H"
#include "BrowningPotter.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{
    defineTypeNameAndDebug(BrowningPotter, 0);
    addToRunTimeSelectionTable
    (
        saturationModel,
        BrowningPotter,
        saturationModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::BrowningPotter::
BrowningPotter
(
    const phaseChangeModel& pcm,
    const dictionary& dict, 
    const objectRegistry& objReg
)
:
    saturationModel
    (
        pcm,
        dict,
        objReg
    ),
    iT_(pcm.pair().iT()),
    p_(pcm.mesh().lookupObject<volScalarField>("p"))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::saturationModels::BrowningPotter::valuePSat
(
    const label& celli
) const
{
    return exp(valueLnPSat(celli));
}

Foam::scalar Foam::saturationModels::BrowningPotter::valuePSatPrime
(
    const label& celli
) const
{
    const scalar& T(iT_[celli]);
    return valuePSat(celli)*(-0.4672/T + 12633.37/sqr(T));
}

Foam::scalar Foam::saturationModels::BrowningPotter::valueLnPSat
(
    const label& celli
) const
{
    const scalar& T(iT_[celli]);
    //- The + log(1e6) is to have p in Pa rather than MPa
    return 
        (11.9463 - 12633.37/T - 0.4672*log(T)) + log(1e6); 
}

Foam::scalar Foam::saturationModels::BrowningPotter::valueTSat
(
    const label& celli
) const
{
    return 
        923840.0/
        (
        -   11275 
        +   Foam::sqrt
            (
                127125625 + 1847680*
                (
                    7.8270 
                -   log(p_[celli]/1e6) // p in MPa
                ) 
            )
        );
}

// ************************************************************************* //
