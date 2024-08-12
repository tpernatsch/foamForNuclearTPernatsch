
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::thermalExpansionSwindemanHastelloyN

Description
    Class modelling thermal expansion of H-N from Swinderman.

SourceFiles
    thermalExpansionSwindemanHastelloyN.C

 \mainauthor
    Alejandra de Lara - University of Cambridge

 \contribution
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour)\n

 \date
    Sep 2023

\*---------------------------------------------------------------------------*/

#ifndef thermalExpansionSwindemanHastelloyN_H
#define thermalExpansionSwindemanHastelloyN_H

#include "thermalExpansionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermalExpansionSwindemanHastelloyN Declaration
\*---------------------------------------------------------------------------*/

class thermalExpansionSwindemanHastelloyN
:
    public thermalExpansionModel
{
    // Private data

        //- Parameters for the Matpro model
        scalar par1;
        scalar par2;
        scalar par3;
       

        //- Perturbation parameter
        scalar perturb;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermalExpansionSwindemanHastelloyN(const thermalExpansionSwindemanHastelloyN&);

        //- Disallow default bitwise assignment
        void operator=(const thermalExpansionSwindemanHastelloyN&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("SwindemanHastelloyN");

    // Constructors

        //- Construct from dictionary
        thermalExpansionSwindemanHastelloyN
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word defaultModel
        );

    //- Destructor
    virtual ~thermalExpansionSwindemanHastelloyN();


    // Member Functions
    
    //- Update thermalExpansion
    virtual void correct(symmTensorField& sf, const scalarField& T, const labelList& addr) const;

    //- Calculate alpha*T
    const symmTensor setAlphaT(const scalar& Ti) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
