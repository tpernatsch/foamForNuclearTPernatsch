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
    Foam::fluidPartitionModels::constant

Description
    Constant interfactial area repartition model

SourceFiles
    constant.C

DictEntryExample ("" is the placeholder for the variable to be inserted)

    fluidPartitionModel
    {
        type        "constant";
        value1      "";
        value2      "";   
    }

    Note that at least one between value1 or value2 is required. If both are
    provided, their sum must be equal to one

\*---------------------------------------------------------------------------*/

#ifndef constant_H
#define constant_H

#include "fluidPartitionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidPartitionModels
{

/*---------------------------------------------------------------------------*\
                           Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public fluidPartitionModel
{
private:

    //- Fraction of interfacial area of the structure that is in contact
    //  with the dispersed phase
    volScalarField fracD_;

    //- Fraction of interfacial area of the structure that is in contact
    //  with the continuous phase
    volScalarField fracC_;

public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct for a fluid-fluid pair
        constant
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fluid& dispersed,
            const fluid& continuous,
            const structureModel& structure
        );

    //- Destructor
    virtual ~constant();


    // Member Functions

        //- Interfacial area
        virtual tmp<volScalarField> fracD() const;

        //- Interfacial area
        virtual tmp<volScalarField> fracC() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidPartitionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //