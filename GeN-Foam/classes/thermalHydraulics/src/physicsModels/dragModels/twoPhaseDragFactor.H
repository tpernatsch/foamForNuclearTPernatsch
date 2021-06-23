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

Description

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

Class
    Foam::twoPhaseDragFactor

SourceFiles
    twoPhaseDragFactor.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseDragFactor_H
#define twoPhaseDragFactor_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FSPair.H"
#include "twoPhaseDragMultiplierModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class twoPhaseDragFactor Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseDragFactor
    : public IOdictionary
{
protected :

    //- Ref to mesh
    const fvMesh& mesh_;

    //- Ref to FSPair
    FSPair& F1SPair_;
    FSPair& F2SPair_;

    //-
    volTensorField KdTotU_;
    
    //- Drag coefficient for an isotropiuc drag factor
    autoPtr<twoPhaseDragMultiplierModel> twoPhaseDragMultiplierPtr_;

public:

    TypeName("twoPhaseDragFactor");

    //- Constructors

        twoPhaseDragFactor
        (
            FSPair& F1SPair,
            FSPair& F2SPair,
            const dictionary& dict
        );

        //- Destructor
        virtual ~twoPhaseDragFactor(){}

    //- Member Functions

        //-
        void correct();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
