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

Class
    Foam::heatTransferModels::nucleateBoiling

Description
    Heat transfer coefficient computed via a super-position model in the form:

    htc = F*htcC + S*htcPB

    where:

    -   htcC is a purely convective heat transfer coefficient
    -   htcPB is a purely pool-boiling related heat transfer coefficient
    -   F is a flow enhancement factor that quantifies the convective HTC
        enhancement due to increased wall turbulence due to vapour bullbe
        departure 
    -   S is a suppression factor that is meant to "turn-off" pool boiling
        as the vapour flow quality increased (as the boiling transitions to
        film boiling).

    In particular, 1 <= F <= inf and 0 <= S <= 1. All these terms are run-time 
    selectable.

Author
    Stefan Radman (stefan.radman@epfl.ch, stefanradman92@gmail.com)

SourceFiles
    nucleateBoiling.C

\*---------------------------------------------------------------------------*/

#ifndef nucleateBoiling_H
#define nucleateBoiling_H

#include "heatTransferModel.H"
#include "flowFactorModel.H"
#include "suppressionFactorModel.H"
#include "poolBoilingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace heatTransferModels
{

using namespace nucleateBoilingSubModels;

/*---------------------------------------------------------------------------*\
                           Class nucleateBoiling Declaration
\*---------------------------------------------------------------------------*/

class nucleateBoiling
:
    public heatTransferModel
{
private:

    mutable autoPtr<flowFactorModel> flowFactor_;

    mutable autoPtr<suppressionFactorModel> suppressionFactor_;

    mutable autoPtr<heatTransferModel> convectionHeatTransfer_;

    mutable autoPtr<poolBoilingModel> poolBoilingHeatTransfer_;

    //-
    mutable volScalarField F_;
    mutable volScalarField S_;

    //- Ref to mass transfer
    const volScalarField dmdt_;

public:

    TypeName("nucleateBoiling");


    //- Constructors

        //- Construct for a fluid-structure pair
        nucleateBoiling
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FSPair& FSPair,
            const wordList& regions
        );

        //- Destructor
        virtual ~nucleateBoiling(){}


    //- Member Functions

        //- Drag coefficient
        virtual void correctHtc(volScalarField& htc) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //