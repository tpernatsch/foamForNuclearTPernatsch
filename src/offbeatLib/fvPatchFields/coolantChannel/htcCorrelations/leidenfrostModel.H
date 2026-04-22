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
    Foam::leidenfrostModel

Description
    Leidenfrost (minimum film boiling) temperature model.

    Models that directly give T_Leidenfrost (empirical):
    - GroeneveldStewartCorrected : Groeneveld-Stewart with quality correction
    - GroeneveldStewart          : Groeneveld-Stewart (1973)
    - Ohnishi                    : Ohnishi et al.
    - Gotovskij                  : Gotovskij

    Models that derive T_Leidenfrost from the minimum film boiling flux:
    - Zuber : Zuber (1958) minimum pool boiling flux.
              T_Leidenfrost is found as T_sat + q_Leid / hFB.
              For a self-consistent result, call the BC's
              solveFilmBoilingTemperature(faceI, qLeid) instead.

SourceFiles
    leidenfrostModel.C

\*---------------------------------------------------------------------------*/

#ifndef leidenfrostModel_H
#define leidenfrostModel_H

#include "scalar.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class leidenfrostModel Declaration
\*---------------------------------------------------------------------------*/

class leidenfrostModel
{
public:

    enum class Model
    {
        GroeneveldStewartCorrected,
        GroeneveldStewart,
        Ohnishi,
        Gotovskij,
        Zuber
    };


    // Constructors

        leidenfrostModel(const word& name);


    // Member Functions

        //- Return the Leidenfrost temperature [K].
        //  hFB is the film boiling HTC evaluated at the current wall
        //  temperature; used only by the Zuber model.
        scalar compute
        (
            scalar T0,  // bulk coolant temperature [K]
            scalar p,   // pressure                 [Pa]
            scalar Q,   // thermodynamic quality    [-]
            scalar hFB  // film boiling HTC         [W/(m²·K)]
        ) const;

        //- True if this model derives T_Leid from a minimum heat flux (Zuber).
        //  When true, call computeMinFlux() then solveFilmBoilingTemperature()
        //  for a self-consistent result.
        bool needsFilmBoilingTemperatureSolve() const;

        //- Return the minimum film boiling flux q_Leid [W/m²] (Zuber only).
        //  Only meaningful when needsFilmBoilingTemperatureSolve() is true.
        scalar computeMinFlux(scalar p) const;

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
