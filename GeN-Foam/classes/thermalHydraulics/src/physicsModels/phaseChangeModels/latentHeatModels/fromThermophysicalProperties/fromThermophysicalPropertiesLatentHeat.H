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

Class
    Foam::latentHeatModels::fromThermophysicalPropertiesLatentHeat

Description
    Class that computes latent heat based on the specified enthalpies of
    formation of the fluids. These are specified under the Hf keyword
    in the thermophysical properties of the two fluids. The latent heat is thus
    computed as Hf.vapour - Hf.liquid.

SourceFiles
    fromThermophysicalPropertiesLatentHeat.C

\*---------------------------------------------------------------------------*/

#ifndef fromThermophysicalPropertiesLatentHeat_H
#define fromThermophysicalPropertiesLatentHeat_H

#include "latentHeatModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace latentHeatModels
{

/*---------------------------------------------------------------------------*\
            Class fromThermophysicalPropertiesLatentHeat Declaration
\*---------------------------------------------------------------------------*/

class fromThermophysicalPropertiesLatentHeat
:
    public latentHeatModel
{
public:

    TypeName("fromThermophysicalProperties");

    //- Constructors

        //- Construct from a dictionary
        fromThermophysicalPropertiesLatentHeat
        (
            const phaseChangeModel& pcm,
            const dictionary& dict, 
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~fromThermophysicalPropertiesLatentHeat(){}

    //- Member Functions

        //-
        virtual void correctField(volScalarField& L) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace latentHeatModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
