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
    Foam::subcooledBoilingHTCModel

Description
    Subcooled nucleate boiling heat transfer correlations.

    Supported models:
    - Chen       : Chen (1966) superposition model
    - JensLottes : Jens-Lottes (1951) empirical correlation
    - Thom       : Thom et al. (1965) empirical correlation
    - Rohsenow   : Rohsenow (1952) pool-boiling correlation

SourceFiles
    subcooledBoilingHTCModel.C

\*---------------------------------------------------------------------------*/

#ifndef subcooledBoilingHTCModel_H
#define subcooledBoilingHTCModel_H

#include "scalar.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class subcooledBoilingHTCModel Declaration
\*---------------------------------------------------------------------------*/

class subcooledBoilingHTCModel
{
public:

    enum class Model
    {
        Chen,
        JensLottes,
        Thom,
        Rohsenow
    };


    // Constructors

        subcooledBoilingHTCModel(const word& name);


    // Member Functions

        //- Compute HTC [W/(m²·K)].
        //  Returns max(hSCB, hSP) to guarantee the boiling term adds to
        //  the single-phase contribution.
        scalar compute
        (
            scalar T0,    // bulk coolant temperature [K]
            scalar Tw,    // wall temperature         [K]
            scalar p,     // pressure                 [Pa]
            scalar G,     // mass flux                [kg/(m²·s)]
            scalar DH,    // hydraulic diameter       [m]
            scalar qWall, // wall heat flux magnitude [W/m²]
            scalar hSP    // single-phase HTC         [W/(m²·K)]
        ) const;

        Model model() const { return model_; }


private:

    Model model_;

    static Model parse(const word& name);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
