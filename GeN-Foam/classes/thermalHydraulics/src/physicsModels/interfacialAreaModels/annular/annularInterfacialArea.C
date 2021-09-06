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
#include "annularInterfacialArea.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interfacialAreaModels
{
    defineTypeNameAndDebug(annular, 0);
    addToRunTimeSelectionTable
    (
        interfacialAreaModel,
        annular,
        interfacialAreaModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfacialAreaModels::annular::annular
(
    const FFPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    interfacialAreaModel
    (
        pair,
        dict,
        objReg
    ),
    dispersed_(pair.alphaDispersed()),
    continuous_(pair.alphaContinuous()),
    DhDispersed_(pair.DhDispersed()),
    cutoffAlpha_(dict.lookupOrDefault<scalar>("cutoffAlpha", 0.9))
{
    if 
    (
        !(pair.fluid1().isLiquid() and pair.fluid2().isGas())
    and !(pair.fluid2().isLiquid() and pair.fluid1().isGas())
    )
    {
        FatalErrorInFunction
            << "Either phase " << pair.fluid1().name() << " or " 
            << pair.fluid2().name()
            << " have an undetermined stateOfMatter (should be specified in "
            << "phaseProperties." << pair.fluid1().name() 
            << "Properties and/or "
            << "phaseProperties." << pair.fluid2().name() << "Properties)"
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::interfacialAreaModels::annular::value
(
    const label& celli
) const
{
    scalar aSum(dispersed_[celli]+continuous_[celli]);
    scalar a(dispersed_[celli]/aSum);
    return 
        aSum*
        (
            (a < cutoffAlpha_) ? 
            4.0*a/DhDispersed_[celli] : 
            4.0*a/DhDispersed_[celli]*(1.0-sqrt(a))/(1.0-cutoffAlpha_)
        );

}

// ************************************************************************* //
