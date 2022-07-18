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
#include "GroeneveldStewartTLF.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LeidenFrostTemperatureModels
{
    defineTypeNameAndDebug(GroeneveldStewart, 0);
    addToRunTimeSelectionTable
    (
        TLFModel,
        GroeneveldStewart,
        LeidenFrostTemperatureModels
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LeidenFrostTemperatureModels::GroeneveldStewart::GroeneveldStewart
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    TLFModel
    (
        pair,
        dict,
        objReg
    ),
    criticalPressure_(dict.get<scalar>("criticalPressure"))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::LeidenFrostTemperatureModels::GroeneveldStewart::value
(
    const label& celli
) const
{
    const scalar& pli(p_[celli]);
    const scalar& Tsati(Tsat_[celli]);

    scalar Tmin(0);
    if (pli<9*1e6) // Correlation from GroeneveldStewart, valid for pressure P<9 MPa
    {
        Tmin = 557.85+44.1*pli*(1e-6)-3.72*pow(pli*1e-6,2);
    }
    else // Ramp up to critical pressure 
    {
        scalar DeltaTmin(557.85+44.1*9-3.72*pow(9,2)-Tsati);
        Tmin = Tsati+(criticalPressure_-pli)/(criticalPressure_-9*1e6)*DeltaTmin;  
    }

    return Tmin;

}

// ************************************************************************* //
