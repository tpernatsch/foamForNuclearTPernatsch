/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
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
    along with OpenFOAM. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::filmBoilingHTCModel

Description
    Post-CHF (film boiling) heat transfer correlations.

    Supported models:
    - GroeneveldAnnular  : Groenveld (1973) annular geometry
    - GroeneveldTubular  : Groenveld (1973) tubular geometry
    - BishopSandbergTong : Bishop-Sandberg-Tong (1965)
    - Frederking         : Frederking (1966) pool film boiling
    - Sakurai            : Sakurai (1990) combined convective/film boiling

SourceFiles
    filmBoilingHTCModel.C

\*---------------------------------------------------------------------------*/

#ifndef filmBoilingHTCModel_H
#define filmBoilingHTCModel_H

#include "scalar.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class filmBoilingHTCModel Declaration
\*---------------------------------------------------------------------------*/

class filmBoilingHTCModel
{
public:

    enum class Model
    {
        GroeneveldAnnular,
        GroeneveldTubular,
        BishopSandbergTong,
        Frederking,
        Sakurai
    };


    // Constructors

        filmBoilingHTCModel(const word& name);


    // Member Functions

        //- Compute film boiling HTC [W/(m²·K)].
        scalar compute
        (
            scalar T0,  // bulk coolant temperature [K]
            scalar Tw,  // wall temperature         [K]
            scalar p,   // pressure                 [Pa]
            scalar Q,   // thermodynamic quality    [-]
            scalar G,   // mass flux                [kg/(m²·s)]
            scalar DH   // hydraulic diameter       [m]
        ) const;

        Model model() const { return model_; }


private:

    Model model_;

    static const scalar g_;

    static Model parse(const word& name);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
