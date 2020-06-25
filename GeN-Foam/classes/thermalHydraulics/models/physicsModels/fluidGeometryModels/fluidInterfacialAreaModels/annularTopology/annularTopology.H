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
    Foam::fluidInterfacialAreaModels::annularTopology

Description
    Interfacial area model that assumes that the dispersed phase flow is
    annular (i.e. cylindrical topology)

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    annularTopology.C

\*---------------------------------------------------------------------------*/

#ifndef annularTopology_H
#define annularTopology_H

#include "fluidInterfacialAreaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidInterfacialAreaModels
{

/*---------------------------------------------------------------------------*\
                           Class annularTopology Declaration
\*---------------------------------------------------------------------------*/

class annularTopology
:
    public fluidInterfacialAreaModel
{
private:

    const volScalarField& DhDispersed_;

    scalar cutoffAlpha_;

public:

    //- Runtime type information
    TypeName("annularTopology");


    // Constructors

        //- Construct for a fluid-fluid pair
        annularTopology
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fluid& dispersed,
            const fluid& continuous
        );

    //- Destructor
    virtual ~annularTopology();


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