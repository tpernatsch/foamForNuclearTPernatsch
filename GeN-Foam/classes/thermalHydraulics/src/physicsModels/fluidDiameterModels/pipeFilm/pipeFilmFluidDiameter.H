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
    Foam::pipeFilm

Description
    Correlation to be used in annular flow scenarios in which the phase of
    interest to which this model is applied flows in a film along the walls
    of a circular pipe, so that the film thickness is proportional to
    1-sqrt(1-alphaN) with alphaN being the normalized phase fraction of the
    phase of interest

Author
    Stefan Radman (stefan.radman@alumni.epfl.ch; sradman@protonmail.com)

SourceFiles
    pipeFilmFluidDiameter.C

\*---------------------------------------------------------------------------*/

#ifndef pipeFilm_H
#define pipeFilm_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidDiameterModel.H"

namespace Foam
{
namespace fluidDiameterModels
{

/*---------------------------------------------------------------------------*\
            Class pipeFilm Declaration
\*---------------------------------------------------------------------------*/

class pipeFilm
:
    public fluidDiameterModel
{
protected:

    //- Minimum allowable value of the phase fraction to avoid null fluid 
    //  diameter
    scalar residualAlpha_;

public:

    //- Runtime type information
    TypeName("pipeFilm");

    // Constructors

        pipeFilm
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~pipeFilm(){}

    // Member Functions

        //- Return fraction of structure interfacial contact in
        //  contact with the dispersed fluid
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidDiameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
