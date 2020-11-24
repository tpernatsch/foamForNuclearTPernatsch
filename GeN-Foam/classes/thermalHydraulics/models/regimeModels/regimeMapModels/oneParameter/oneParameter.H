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
    Foam::regimeMapModel::oneParameter

Description
    Regime map model that depend on a single scalar parameter value.

Autor
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    oneParameter.C

\*---------------------------------------------------------------------------*/

#ifndef oneParameter_H
#define oneParameter_H

#include "regimeMapModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace regimeMapModels
{

/*---------------------------------------------------------------------------*\
                           Class oneParameter Declaration
\*---------------------------------------------------------------------------*/

class oneParameter
:
    public regimeMapModel
{
public:

    enum class interpolationMode
    {
        linear, 
        quadratic
    };

    static const Enum<interpolationMode> interpolationModeNames_;

private:

    const scalarField& parameter_;

    wordList orderedRegimeNames_;

    scalarList thresholds_;

    const interpolationMode interpolationMode_;

public:

    //- Runtime type information
    TypeName("oneParameter");

    //- Constructors

        oneParameter
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& physicsModelsDict
        );

        //- Destructor
        virtual ~oneParameter(){}

    //- Member Functions

        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regimeMapModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //