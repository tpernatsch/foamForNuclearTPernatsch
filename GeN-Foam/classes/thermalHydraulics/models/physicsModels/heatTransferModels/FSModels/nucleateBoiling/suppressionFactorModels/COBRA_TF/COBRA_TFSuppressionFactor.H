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
    Suppression factor used in the COBRA-TF code, in which 
    S = (Tfluid-Tsat)/(Twall-Tsat)

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    COBRA_TFSuppressionFactor.C

\*---------------------------------------------------------------------------*/

#ifndef COBRA_TF_H
#define COBRA_TF_H

#include "suppressionFactorModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declaration
class structureModel;

namespace heatTransferModels
{

namespace nucleateBoilingSubModels
{

namespace suppressionFactorModels
{

/*---------------------------------------------------------------------------*\
                           Class COBRA_TF Declaration
\*---------------------------------------------------------------------------*/

class COBRA_TF
:
    public suppressionFactorModel
{
protected :

    //- Ref to structure wall temperature
    const volScalarField& Twall_;

    //- Ref to fluid temperature
    const volScalarField& Tf_;

    //- Ref to saturation temperature
    const volScalarField& Tsat_;


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

        //- Correct the suppression factor
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace suppressionFactorModels
} // End namespace nucleateBoilingSubModels
} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //