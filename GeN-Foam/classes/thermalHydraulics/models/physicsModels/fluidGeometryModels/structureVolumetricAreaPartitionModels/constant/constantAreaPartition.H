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
    Foam::structureVolumetricAreaPartitionModels::constant

Description
    Constant interfactial area repartition model

SourceFiles
    constant.C

DictEntryExample ("" is the placeholder for the variable to be inserted)

    structureVolumetricAreaPartitionModel
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

#include "structureVolumetricAreaPartitionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace structureVolumetricAreaPartitionModels
{

/*---------------------------------------------------------------------------*\
                           Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public structureVolumetricAreaPartitionModel
{
private:

    //- Fraction of interfacial area of the structure that is in contact
    //  with fluid1
    volScalarField frac1_;

    //- Fraction of interfacial area of the structure that is in contact
    //  with fluid2
    volScalarField frac2_;

public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct for a fluid-fluid pair
        constant
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fluid& fluid1,
            const fluid& fluid2,
            const structureModel& structure
        );

    //- Destructor
    virtual ~constant();


    // Member Functions

        //- Interfacial area
        virtual tmp<volScalarField> frac1() const;

        //- Interfacial area
        virtual tmp<volScalarField> frac2() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace structureVolumetricAreaPartitionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //