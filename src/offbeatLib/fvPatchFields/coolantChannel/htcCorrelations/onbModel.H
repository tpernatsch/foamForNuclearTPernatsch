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
    Foam::onbModel

Description
    Onset of Nucleate Boiling (ONB) temperature model.

    Returns the wall temperature at which nucleate boiling is initiated.

    Supported models:
    - Basu           : Basu (2002) – physics-based ONB criterion
    - BerglesRohsenow: Bergles and Rohsenow (1964) – iterative empirical criterion

SourceFiles
    onbModel.C

\*---------------------------------------------------------------------------*/

#ifndef onbModel_H
#define onbModel_H

#include "scalar.H"
#include "word.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class onbModel Declaration
\*---------------------------------------------------------------------------*/

class onbModel
{
public:

    enum class Model
    {
        Basu,
        BerglesRohsenow
    };


    // Constructors

        onbModel(const word& name);


    // Member Functions

        //- Return ONB temperature [K].
        //  hSP is the local single-phase HTC [W/(m²·K)].
        scalar compute
        (
            scalar T0,        // bulk coolant temperature              [K]
            scalar p,         // pressure                              [Pa]
            scalar kCoolant,  // coolant thermal conductivity at film T [W/(m·K)]
            scalar hSP,       // single-phase HTC                      [W/(m²·K)]
            label  faceI      // face index (for warning messages only)
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
