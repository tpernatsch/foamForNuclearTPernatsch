
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
    Foam::YoungModulusModel

Description
    Mother class for Young's Modulus models.

SourceFiles
    YoungModulusModel.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef YoungModulusModel_H
#define YoungModulusModel_H

#include "primitiveFields.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class YoungModulusModel Declaration
\*---------------------------------------------------------------------------*/

class YoungModulusModel
{
    // Private data
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        YoungModulusModel(const YoungModulusModel&);

        //- Disallow default bitwise assignment
        void operator=(const YoungModulusModel&);

protected:
    
    // Protected data
        
        //- Reference to object mesh
        const fvMesh& mesh_;

public:

    //- Runtime type information
    TypeName("YoungModulusModel");

    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        YoungModulusModel,
        dictionary,
        (
            const fvMesh& mesh, 
            const dictionary& dict, 
            const word defaultModel
        ),
        (mesh, dict, defaultModel)
    );

    // Constructors

        //- Construct from dictionary
        YoungModulusModel
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel   
        );


    // Selectors

        //- Select from dictionary
        static autoPtr<YoungModulusModel> New
        (   
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel 
        );

    //- Destructor
    virtual ~YoungModulusModel();


    // Member Functions
    
    //- Update Young Modulus  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
