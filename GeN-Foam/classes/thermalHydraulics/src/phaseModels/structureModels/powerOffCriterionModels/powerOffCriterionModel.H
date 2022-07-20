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
    Foam::powerOffCriterionModel

Description


Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    powerOffCriterionModel.C
    newFluidDiameterModel.C

\*---------------------------------------------------------------------------*/

#ifndef powerOffCriterionModel_H
#define powerOffCriterionModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class powerOffCriterionModel Declaration
\*---------------------------------------------------------------------------*/

class powerOffCriterionModel
:
    public IOdictionary
{
protected:

    //- Ref to fluid
    const fvMesh& mesh_;

public:

    //- Runtime type information
    TypeName("powerOffCriterionModel");

    //- Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            powerOffCriterionModel,
            powerOffCriterionModels,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );

    //- Constructors

        powerOffCriterionModel
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Destructor
        virtual ~powerOffCriterionModel(){}

    //- Selectors

        static autoPtr<powerOffCriterionModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Member Functions

        //- Return true if power should be turned off
        virtual bool powerOffCriterion() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
