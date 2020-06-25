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
    Foam::powerOffControlModels::fieldValue

Description
    Turns power off if the field value maximum exceeds a provided value.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    fieldValue.C

DictEntryExample ("" is the placeholder for the variable to be inserted)

    powerOffControlModel
    {
        type        "fieldValue";
        value       "";   
    }


\*---------------------------------------------------------------------------*/

#ifndef fieldValue_H
#define fieldValue_H

#include "powerOffCriterionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace powerOffCriterionModels
{

/*---------------------------------------------------------------------------*\
                           Class fieldValue Declaration
\*---------------------------------------------------------------------------*/

class fieldValue
:
    public powerOffCriterionModel
{
public:

    enum fieldOp
    {
        max,
        min
    };

    static const Enum<fieldOp> fieldOpNames_;

    enum criterion
    {
        above,
        below
    };

    static const Enum<criterion> criterionNames_;

private:

    const word fieldName_;

    const fieldOp fieldOp_;

    const criterion criterion_;

    const scalar threshold_;

    const scalar timeDelay_;

    scalar time0_;

    const volScalarField* fieldPtr_;

    bool InfoFlag_;

public:

    //- Runtime type information
    TypeName("fieldValue");

    //- Constructors

        fieldValue
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~fieldValue();

    //- Member Functions

        //- Interfacial area
        virtual bool powerOffCriterion();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace powerOffControlModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //