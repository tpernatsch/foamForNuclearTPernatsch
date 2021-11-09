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
    Foam::latentHeatModels::water

Description
    This class describes the latent heat of vaporization of water in the 
    0.01-350 C range. It was obtained from a fit (in the form of 
    L = A+B*log(C-T) with T being the liquid temperature (K)) of NIST water 
    data at saturation found at:
    https://www.nist.gov/system/files/documents/srd/NISTIR5078-Tab1.pdf

    It performs reasonably well up to ~ 365 C as well, considering that the
    critical water temperature is 375.16 C (at which point the latent heat 
    becomes 0, but is limited for obvious numerical reasons as this code is not 
    meant for sub-supercritical fluid transition simulations)

Author
    Stefan Radman (sradman@pm.me / stefan.radman@epfl.ch)

SourceFiles
    waterLatentHeat.C

\*---------------------------------------------------------------------------*/

#ifndef waterLatentHeat_H
#define waterLatentHeat_H

#include "latentHeatModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace latentHeatModels
{

/*---------------------------------------------------------------------------*\
            Class water Declaration
\*---------------------------------------------------------------------------*/

class water
:
    public latentHeatModel
{
public:

    TypeName("water");

    //- Constructors

        //- Construct from a dictionary
        water
        (
            const phaseChangeModel& pcm,
            const dictionary& dict, 
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~water(){}

    //- Member Functions

        virtual scalar value(const label& celli) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace latentHeatModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
