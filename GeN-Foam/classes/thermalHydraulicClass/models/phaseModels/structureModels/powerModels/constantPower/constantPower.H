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
    Foam::structureModels::constantPower

Description
    Describes a generic structure with a constant internal power production.
    Solves an energy equation for the structure, updating the average
    structure temperature (assumed equal to the structure surface temperature
    for fluid heat transfer purposes).

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    constantPower.C

\*---------------------------------------------------------------------------*/

#ifndef constantPower_H
#define constantPower_H

#include "powerModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declaration
class structureModel;

namespace powerModels
{

/*---------------------------------------------------------------------------*\
                           Class constantPower Declaration
\*---------------------------------------------------------------------------*/

class constantPower
:
    public powerModel
{
protected :

    //- Volumetric power
    volScalarField powerDensity_;

    //- Structure temperature
    volScalarField T_;

    //- Heat capacity of the structure
    dimensionedScalar Cp_;

    //- Structure density
    dimensionedScalar rho_;

public:

    //- Runtime type information
    TypeName("constantPower");

    //- Constructors

        constantPower
        (
            const structureModel& structure,
            const dictionary& dict,
            const word name
        );

        //- Destructor
        virtual ~constantPower();

    //- Member Functions

        //- Return surface temperature
        virtual tmp<volScalarField> T() const;

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
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace powerModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //