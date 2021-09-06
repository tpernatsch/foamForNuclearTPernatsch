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
    Foam::LockhartMartinelli

Description
    Class to handle fluid-structure friction coefficient (drag) multipliers.

Author
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    LockhartMartinelli.C
    newLockhartMartinelli.C

\*---------------------------------------------------------------------------*/

#ifndef LockhartMartinelli_H
#define LockhartMartinelli_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "twoPhaseDragMultiplierModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace twoPhaseDragMultiplierModels
{

/*---------------------------------------------------------------------------*\
                    Class LockhartMartinelli Declaration
\*---------------------------------------------------------------------------*/

class LockhartMartinelli
:
    public twoPhaseDragMultiplierModel
{
private:

    //- C in phi2 = 1 + C/X + 1/X^2
    scalar C_;

public:

    //- Runtime type information
    TypeName("LockhartMartinelli");

    // Constructors

        LockhartMartinelli
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~LockhartMartinelli(){}

    // Member Functions

        //- Return multiplier
        virtual scalar phi2(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoPhaseDragMultiplierModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
