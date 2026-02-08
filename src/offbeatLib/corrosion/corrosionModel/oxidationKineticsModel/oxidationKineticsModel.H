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
    Foam::oxidationKineticsModel

Description
    Base class for oxidation kinetics model

SourceFiles
    oxidationKineticsModel.C

\*---------------------------------------------------------------------------*/

#ifndef oxidationKineticsModel_H
#define oxidationKineticsModel_H

#include "runTimeSelectionTables.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class oxidationKineticsModel Declaration
\*---------------------------------------------------------------------------*/

class oxidationKineticsModel
{

public:

    //- Runtime type information
    TypeName("oxidationKineticsModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            oxidationKineticsModel,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );

    // Constructors

        //- Construct from dictionary
        oxidationKineticsModel
        (
            const dictionary& dict
        );

    // Selectors

        //- Select from dictionary
        static autoPtr<oxidationKineticsModel> New
        (
            const word& type,
            const dictionary& dict
        );

    //- Destructor
    virtual ~oxidationKineticsModel();


    // Member Functions

        //- Correct the oxide layer thickness [m]
        virtual void correctOxideThickness
        (
            scalarField& S,
            const scalarField& S0,
            const scalarField& Ti,
            const scalarField& Tb,
            const scalarField& phi,
            const scalar dt
        ) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
