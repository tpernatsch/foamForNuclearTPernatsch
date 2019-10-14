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
    Foam::structureModels::constantFluidPower

Description
    Constant power directly added to the fluid in the region where the
    structure exists. Obviously unphysical, for development purposes only.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    constantFluidPower.C

\*---------------------------------------------------------------------------*/

#ifndef constantFluidPower_H
#define constantFluidPower_H

#include "powerModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace powerModels
{

/*---------------------------------------------------------------------------*\
                           Class constantFluidPower Declaration
\*---------------------------------------------------------------------------*/

class constantFluidPower
:
    public powerModel
{
protected :

    //- Volumetric power
    volScalarField& P_;

public:

    //- Runtime type information
    TypeName("constantFluidPower");

    //- Constructors

        constantFluidPower
        (
            const structureModel& structure,
            const dictionary& dict,
            const word name,
            volScalarFieldTable& fields
        );

        //- Destructor
        virtual ~constantFluidPower();

    //- Member Functions

        //- In two-phase solvers :
        //  HT = H1*frac1*T1 + H2*frac2*T2
        //  H =  H1*frac1 + H2*frac2
        //- In a mono-phase solver :
        //  HT = H*T
        //  H = H
        virtual void correct
        (
            const volScalarField& HT,
            const volScalarField& H
        ) {}

        //- Return fluid heat source given fluid temperature and fluid heat
        //  transfer coefficient
        virtual tmp<volScalarField> heatSource
        (
            const volScalarField& Tfluid,
            const volScalarField& h
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace powerModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //