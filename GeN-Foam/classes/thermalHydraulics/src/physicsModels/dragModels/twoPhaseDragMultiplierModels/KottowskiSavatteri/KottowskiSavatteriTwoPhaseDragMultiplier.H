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
    Foam::KottowskiSavatteri

Description
    Model by Kottowski and Savatteri for the multiplier, specifically for
    sodium two-phase flows:

    phi2 = 10^(2*( 0.1046*(log10X)^2 - 0.5098*log10X + 0.6252 ))

    log10 is the log in basis 10, while X2 is the well-known 
    Lockhart-Martinelli parameter and is defined as:

    X2 = (mu_I/mu_O)^0.2 * (flowQuality_I/flowQuality_O)^1.8 * rho_O/rho_I

    Please note that X in use in phi2 is X = sqrt(X2). The correlation for
    phi2 was derived in the data range 7e-2 < X < 30.

    For futher info see H.M. Kottowsi, C. Savatteri, "Fundamentals of Liquid
    Metal Boiling Thermohydraulics", Nucl. Eng. Des. vol. 82, pp. 281-304, 1984

Author
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    KottowskiSavatteri.C

\*---------------------------------------------------------------------------*/

#ifndef KottowskiSavatteriTwoPhaseDragMultiplier_H
#define KottowskiSavatteriTwoPhaseDragMultiplier_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "twoPhaseDragMultiplierModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace twoPhaseDragMultiplierModels
{

/*---------------------------------------------------------------------------*\
                    Class KottowskiSavatteri Declaration
\*---------------------------------------------------------------------------*/

class KottowskiSavatteri
:
    public twoPhaseDragMultiplierModel
{
public:

    //- Runtime type information
    TypeName("KottowskiSavatteri");

    // Constructors

        KottowskiSavatteri
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~KottowskiSavatteri(){}

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
