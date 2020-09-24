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
    Foam::structureModels::COBRA_TF

Description
    Describes a generic structure with a constant internal power production.
    Solves an energy equation for the structure, updating the average
    structure temperature (assumed equal to the structure surface temperature
    for fluid heat transfer purposes).

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    COBRA_TFFlowFactor.C

\*---------------------------------------------------------------------------*/

#ifndef COBRA_TF_H
#define COBRA_TF_H

#include "flowFactorModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declaration
class structureModel;

namespace heatTransferModels
{

namespace nucleateBoilingSubModels
{

namespace flowFactorModels
{

/*---------------------------------------------------------------------------*\
                           Class COBRA_TF Declaration
\*---------------------------------------------------------------------------*/

class COBRA_TF
:
    public flowFactorModel
{
protected :

    //- Superficial fluid Reynolds
    const volScalarField& Re1p_;

    //- Two-phase Reynolds
    const volScalarField& Re2p_;


public:

    //- Runtime type information
    TypeName("COBRA-TF");

    //- Constructors

        COBRA_TF
        (
            const heatTransferModel& htm,
            const objectRegistry& objReg,
            const FSPair& FSPair
        );

        //- Destructor
        virtual ~COBRA_TF(){}

    //- Member Functions

        //- Correct the flow factor
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace flowFactorModels
} // End namespace nucleateBoilingSubModels
} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //