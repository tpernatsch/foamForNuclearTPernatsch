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
    Foam::Kaiser74

Description
    Model by Kaiser et al. for the multiplier, specifically for sodium
    two-phase flows:

    phi = 8.2/(X^0.55)

    X2 is the well-known Lockhart-Martinelli
    parameter defined as:

    X2 = (mu_I/mu_O)^0.2 * (flowQuality_I/flowQuality_O)^1.8 * rho_O/rho_I

    Please note that X in use in phi is X = sqrt(X2) and that phi2 = phi^2

    Kaiser, A., Peppler, W., Vöröss, L., 1974. "Type of flow, pressure drop, 
    and critical heat flux of a two-phase sodium flow". Nucl. Eng. Des. 30, 
    pp. 305–315.

Author
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    Kaiser74.C

\*---------------------------------------------------------------------------*/

#ifndef Kaiser74TwoPhaseDragMultiplier_H
#define Kaiser74TwoPhaseDragMultiplier_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "twoPhaseDragMultiplierModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace twoPhaseDragMultiplierModels
{

/*---------------------------------------------------------------------------*\
                    Class Kaiser74 Declaration
\*---------------------------------------------------------------------------*/

class Kaiser74
:
    public twoPhaseDragMultiplierModel
{
public:

    //- Runtime type information
    TypeName("Kaiser74");

    // Constructors

        Kaiser74
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~Kaiser74(){}

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
