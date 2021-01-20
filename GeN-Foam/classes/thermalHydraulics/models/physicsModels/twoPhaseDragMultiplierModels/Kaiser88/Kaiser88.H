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
    Foam::twoPhaseDragMultiplierModels::Kaiser88

Description
    Model by Kaiser et al. for the multiplier, specifically for sodium
    two-phase flows:

    phi2 = exp(2*(  1.48 - 1.05*log(sqrt(X)) + 0.09*(log(sqrt(X)))^2  ))

    log is the natural logarithm and X2 is the well-known Lockhart-Martinelli
    parameter defined as:

    X2 = (mu_I/mu_O)^0.2 * (flowQuality_I/flowQuality_O)^1.8 * rho_O/rho_I

    Please note that X in use in phi2 is X = sqrt(X2).

    Kaiser, A., Huber, F., Bottoni, M., Dorr, B., 1988. "Contribution to sodium
    boiling heat transfer, pressure drop, and void distribution in multi-pin 
    geometry". In: Proc. 13th Liquid Metal Boiling Working Group (LMBWG), pp. 
    71–98.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    Kaiser88.C

DictEntryExample

    twoPhaseDragMultiplierModel
    {
        type            Kaiser88;
        multiplierFluid liquid;
    }


\*---------------------------------------------------------------------------*/

#ifndef Kaiser88_H
#define Kaiser88_H

#include "twoPhaseDragMultiplierModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace twoPhaseDragMultiplierModels
{

/*---------------------------------------------------------------------------*\
                           Class Kaiser88 Declaration
\*---------------------------------------------------------------------------*/

class Kaiser88
:
    public twoPhaseDragMultiplierModel
{
private:

public:

    //- Runtime type information
    TypeName("Kaiser88");

    //- Constructors

        Kaiser88
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~Kaiser88();

    //- Member Functions

        //- Calculate phi2_;
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoPhaseDragMultiplierModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //