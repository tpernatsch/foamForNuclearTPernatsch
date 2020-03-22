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
    Base abstract class for the fluid and structure classes.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

Class
    Foam::phaseBase

SourceFiles
    phaseBase.C

\*---------------------------------------------------------------------------*/


#ifndef phaseBase_H
#define phaseBase_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "typedefs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseBase Declaration
\*---------------------------------------------------------------------------*/

class phaseBase
:
    public volScalarField
{
protected:

    //- Phase dictionary
    const dictionary dict_;

    //- Reference to mesh
    const fvMesh& mesh_;

    //- Name of phase
    const word name_;

    //- Return the residual phase-fraction for given phase
    //  Used to stabilize the phase momentum as the phase-fraction -> 0
    const dimensionedScalar residualAlpha_;
    

public:

    // Constructors

        phaseBase
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& name,
            bool readIfPresFlag = true,
            bool writeFlag = true
        );

        //- Destructor
        virtual ~phaseBase();

    // Member Functions

        //- To make it pure abstract
        virtual void abstract() = 0;

        // Access

            const fvMesh& mesh() const
            {
                return mesh_;
            }
            const word& name() const
            {
                return name_;
            }
            const dictionary& dict() const
            {
                return dict_;
            }
            const dimensionedScalar& residualAlpha() const
            {
                return residualAlpha_;
            } 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
