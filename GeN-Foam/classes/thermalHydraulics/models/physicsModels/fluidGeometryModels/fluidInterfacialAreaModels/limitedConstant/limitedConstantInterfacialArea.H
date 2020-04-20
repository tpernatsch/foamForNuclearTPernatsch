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
    Foam::fluidInterfacialAreaModels::limitedConstant

Description
    Constant interfacial area that has a prescribed value only in 
    two-phase flow domains, 0 elsewhere

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    limitedConstantInterfacialArea.C

DictEntryExample ("" is the placeholder for the variable to be inserted)

    fluidInterfacialAreaModel
    {
        type        "limitedConstant";
        value       ""; // interfacial area value
    }

\*---------------------------------------------------------------------------*/

#ifndef limitedConstant_H
#define limitedConstant_H

#include "fluidInterfacialAreaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidInterfacialAreaModels
{

/*---------------------------------------------------------------------------*\
                           Class limitedConstant Declaration
\*---------------------------------------------------------------------------*/

class limitedConstant
:
    public fluidInterfacialAreaModel
{
private:

    //- Phase fraction of vapour phase
    const volScalarField iA_;

    //- Residual value of phaseFraction below which iA_ is set to 0;
    const dimensionedScalar residualAlpha_;

public:

    //- Runtime type information
    TypeName("limitedConstant");


    // Constructors

        //- Construct for a fluid-fluid pair
        limitedConstant
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fluid& dispersed,
            const fluid& continuous
        );

    //- Destructor
    virtual ~limitedConstant();


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