
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::conductivitySwindemanHastelloyN

Description
    Model for conductivity of H-N, derived from the work of Swinderman.
    
SourceFiles
    conductivitySwindemanHastelloyN.C

 \mainauthor
    Alejandra de Lara - University of Cambridge

 \contribution
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour)\n

 \date
     Sep 2023

\*---------------------------------------------------------------------------*/

#ifndef conductivitySwindemanHastelloyN_H
#define conductivitySwindemanHastelloyN_H

#include "conductivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class conductivitySwindemanHastelloyN Declaration
\*---------------------------------------------------------------------------*/

class conductivitySwindemanHastelloyN
:
    public conductivityModel
{
    // Private data

        //- Parameters for the  model
        scalar par1;
        scalar par2;

        //- Perturbation factor
        scalar perturb;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        conductivitySwindemanHastelloyN(const conductivitySwindemanHastelloyN&);

        //- Disallow default bitwise assignment
        void operator=(const conductivitySwindemanHastelloyN&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("SwindemanHastelloyN");

    // Constructors

        //- Construct from dictionary
        conductivitySwindemanHastelloyN
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word defaultModel
        );

    //- Destructor
    virtual ~conductivitySwindemanHastelloyN();


    // Member Functions
    
    //- Update conductivity
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
