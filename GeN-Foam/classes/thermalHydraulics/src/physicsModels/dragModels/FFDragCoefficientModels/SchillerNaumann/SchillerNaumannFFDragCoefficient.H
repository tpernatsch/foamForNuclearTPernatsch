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
    Foam::FFDragCoefficientModels::SchillerNaumann

Description
    Drag coefficient in the form fd = coeff*Re^exp with 
    Kd = 0.5*fd*alpha*rho*magU/Dh

Dictionary entries (mandatory)
    coeff   0.316;
    exp     -0.25;

SourceFiles
    SchillerNaumann.C

\*---------------------------------------------------------------------------*/

#ifndef SchillerNaumannFFDragCoefficient_H
#define SchillerNaumannFFDragCoefficient_H

#include "FFDragCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class FFPair;

namespace FFDragCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class SchillerNaumann Declaration
\*---------------------------------------------------------------------------*/

class SchillerNaumann
:
    public FFDragCoefficientModel
{
public:

    TypeName("SchillerNaumann");

    //- Constructors

        //- Construct for a fluid-structure pair
        SchillerNaumann
        (
            const FFPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~SchillerNaumann(){}


    //- Member Functions

        //- Drag coefficient
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FFDragCoefficientModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //