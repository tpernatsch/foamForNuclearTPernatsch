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
    Foam::fluidInterfacialAreaModels::integralConstant

Description
    integralConstant interfacial area model between the two fluid. Only for
    development/verification purposes

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    integralConstantInterfacialArea.C

DictEntryExample ("" is the placeholder for the variable to be inserted)

    fluidInterfacialAreaModel
    {
        type        "integralConstant";
        value       ""; // interfacial area value
    }

\*---------------------------------------------------------------------------*/

#ifndef integralConstant_H
#define integralConstant_H

#include "fluidInterfacialAreaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidInterfacialAreaModels
{

/*---------------------------------------------------------------------------*\
                           Class integralConstant Declaration
\*---------------------------------------------------------------------------*/

class integralConstant
:
    public fluidInterfacialAreaModel
{
private:

    //- Integral interfacial area (not a density!)
    const dimensionedScalar intA0_;

    //- Residual value of phaseFraction below which iA_ is set to 0;
    const dimensionedScalar residualAlpha_;

public:

    //- Runtime type information
    TypeName("integralConstant");


    // Constructors

        //- Construct for a fluid-fluid pair
        integralConstant
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fluid& dispersed,
            const fluid& continuous
        );

    //- Destructor
    virtual ~integralConstant();


    // Member Functions

        //- Interfacial area
        virtual tmp<volScalarField> iA() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidInterfacialAreaModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //