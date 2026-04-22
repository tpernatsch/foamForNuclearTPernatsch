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
    Foam::transitionBoilingHTCModel

Description
    Transition boiling heat transfer correlations (between CHF and Leidenfrost
    points).

    Supported models:
    - McDonoughMilichKing : McDonough-Milich-King (1961)
    - None                : Falls back to the film boiling HTC (lower bound)

SourceFiles
    transitionBoilingHTCModel.C

\*---------------------------------------------------------------------------*/

#ifndef transitionBoilingHTCModel_H
#define transitionBoilingHTCModel_H

#include "scalar.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class transitionBoilingHTCModel Declaration
\*---------------------------------------------------------------------------*/

class transitionBoilingHTCModel
{
public:

    enum class Model
    {
        McDonoughMilichKing,
        None
    };


    // Constructors

        transitionBoilingHTCModel(const word& name);


    // Member Functions

        //- Compute transition boiling HTC [W/(m²·K)].
        //  Returns max(hTB, hFB) so the result never falls below film boiling.
        scalar compute
        (
            scalar p,     // pressure                [Pa]
            scalar qWall, // wall heat flux magnitude [W/m²]
            scalar TCHF,  // CHF wall temperature    [K]
            scalar qCHF,  // critical heat flux      [W/m²]
            scalar hFB    // film boiling HTC        [W/(m²·K)]
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
