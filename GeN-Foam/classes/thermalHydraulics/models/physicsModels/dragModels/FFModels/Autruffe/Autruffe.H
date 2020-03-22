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
    Foam::dragModels::Autruffe

Description
    Autruffe correlation for the pressure drop. To be used
    only for momentum transfer between a solid phase (i.e. stationary) 
    and a fluid phase (i.e. moving).
    By Stefan Radman (stefan.radman@epfl.ch)

SourceFiles
    Autruffe.C

\*---------------------------------------------------------------------------*/

#ifndef Autruffe_H
#define Autruffe_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class Autruffe Declaration
\*---------------------------------------------------------------------------*/

class Autruffe
:
    public dragModel
{
protected:

    //- The FFPair does not tell me which phase is liquid and which is vapour,
    //  so I need these refs

    //- Ref to vapour
    const fluid& vapour_;

    //- Ref to liquid
    const fluid& liquid_;

public:

    TypeName("Autruffe");

    //- Constructors

        //- Construct from a dictionary and a phase pair
        Autruffe
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FFPair& FFPair
        );

        //- Destructor
        virtual ~Autruffe(){}

    //- Member Functions

        //- Drag coefficient
        virtual void correctKd(volTensorField& Kd) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
