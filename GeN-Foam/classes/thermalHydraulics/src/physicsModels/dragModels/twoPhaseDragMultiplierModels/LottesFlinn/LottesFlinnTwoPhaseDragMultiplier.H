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
    Foam::LottesFlinn

Description
    Model in the spirit of Lottes-Flinn, i.e.: 

    phi2 = 1/(alpha^exp) 
    
    with alpha being the normalized phase fraction (by the void fraction) of
    the multiplierFluid and exp defaulting to 2 (which is the actual
    Lottes-Flinn model). The exp can be user-provided (it is know from e.g.
    TRACE development that this exponent can be slightly smaller for bubbly
    flows, e.g. in the 1.7-1.8 range). 

Author
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    LottesFlinn.C
    newLottesFlinn.C

\*---------------------------------------------------------------------------*/

#ifndef LottesFlinnTwoPhaseDragMultiplier_H
#define LottesFlinnTwoPhaseDragMultiplier_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "twoPhaseDragMultiplierModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace twoPhaseDragMultiplierModels
{

/*---------------------------------------------------------------------------*\
                    Class LottesFlinn Declaration
\*---------------------------------------------------------------------------*/

class LottesFlinn
:
    public twoPhaseDragMultiplierModel
{
private:

    //- exp in phi2 = 1/(alpha^exp)
    scalar exp_;

public:

    //- Runtime type information
    TypeName("LottesFlinn");

    // Constructors

        LottesFlinn
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~LottesFlinn(){}

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
