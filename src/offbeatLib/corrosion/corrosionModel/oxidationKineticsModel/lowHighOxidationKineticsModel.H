/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::lowHighOxidationKineticsModel

Description
    Base class for low and high temperature oxidation kinetics models,
    templated for performance.

SourceFiles
    lowHighOxidationKineticsModel.C

\*---------------------------------------------------------------------------*/

#ifndef lowHighOxidationKineticsModel_H
#define lowHighOxidationKineticsModel_H

#include "oxidationKineticsModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class lowHighOxidationKineticsModel Declaration
\*---------------------------------------------------------------------------*/

template<class LowTemperatureModel,class HighTemperatureModel>
class lowHighOxidationKineticsModel
:
    public oxidationKineticsModel
{

    // Private data

        //- Low Temperature oxidation model
        LowTemperatureModel lowTModel_;

        //- High Temperature oxidation model
        HighTemperatureModel highTModel_;

public:

    //- Runtime type information
    TypeName("lowHighOxidationKineticsModel");

    // Constructors

        //- Construct from dictionary
        lowHighOxidationKineticsModel
        (
            const dictionary& dict
        );

    //- Destructor
    ~lowHighOxidationKineticsModel();


    // Member Functions

        //- Correct the oxide layer thickness [m]
        void correctOxideThickness
        (
            scalarField& S,
            const scalarField& S0,
            const scalarField& Ti,
            const scalarField& Tb,
            const scalarField& phi,
            const scalar dt
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "lowHighOxidationKineticsModel.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
