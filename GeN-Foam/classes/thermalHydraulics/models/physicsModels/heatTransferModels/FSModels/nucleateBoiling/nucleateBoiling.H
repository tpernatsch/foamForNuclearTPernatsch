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
    Heat transfer coefficient computed from Nusselt via a correlation in the
    form: A + B*Reynolds^(C)*Prandtl^(D).

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

    //- Ref to structure wall temperature
    const volScalarField& Twall_;

    //- Ref to fluid temperature
    const volScalarField& Tf_;

    //- Ref to saturation temperature
    const volScalarField& Tsat_;

    //- Ref to pressure
    const volScalarField& p_;

    //- Critical pressure in Pa
    scalar pCrit_;

    mutable autoPtr<flowFactorModel> flowFactor_;

    mutable autoPtr<suppressionFactorModel> suppressionFactor_;

    mutable autoPtr<heatTransferModel> convectionHeatTransfer_;

    mutable autoPtr<poolBoilingModel> poolBoilingHeatTransfer_;

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