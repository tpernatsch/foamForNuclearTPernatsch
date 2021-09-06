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
    Foam::FFDragCoefficientModels::Autruffe

Description
    Autruffe correlation for the friction factor between two fluids,
    originally developed for sodium liquid-vapour pairs. See
    Autruffe, M.I., Wilson, G.J., Stewart, B., Kazimi, M.S., 1979. A proposed 
    momentum exchange coefficient for two-phase modeling of sodium boiling. 
    In: Proc. Int. Meeting on Fast Reactor Safety Tech. 4, Seattle, WA, pp. 
    2515–2521

Author
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    Autruffe.C

\*---------------------------------------------------------------------------*/

#ifndef AutruffeFFDragCoefficient_H
#define AutruffeFFDragCoefficient_H

#include "FFDragCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class FFPair;

namespace FFDragCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class Autruffe Declaration
\*---------------------------------------------------------------------------*/

class Autruffe
:
    public FFDragCoefficientModel
{
protected:

    const fluid& vapour_;

public:

    TypeName("Autruffe");

    //- Constructors

        //- Construct for a fluid-structure pair
        Autruffe
        (
            const FFPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~Autruffe(){}


    //- Member Functions

        //- Drag coefficient
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FFDragCoefficientModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //