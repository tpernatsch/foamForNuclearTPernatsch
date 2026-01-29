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
    Foam::YoungModulusHofmanD9

Description
    Class modelling the Young's modulus of D9 Steel cladding using the Hofman 
    correlation.

SourceFiles
    YoungModulusHofmanD9.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef YoungModulusHofmanD9_H
#define YoungModulusHofmanD9_H

#include "YoungModulusModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class YoungModulusHofmanD9 Declaration
\*---------------------------------------------------------------------------*/

class YoungModulusHofmanD9
:
    public YoungModulusModel
{
    // Private data

        //- Parameters for the Matpro model
        scalar par1;
        scalar par2;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        YoungModulusHofmanD9(const YoungModulusHofmanD9&);

        //- Disallow default bitwise assignment
        void operator=(const YoungModulusHofmanD9&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("SteelD9Hofman");

    // Constructors

        //- Construct from dictionary
        YoungModulusHofmanD9
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict, 
            const word defaultModel   
        );
    //- Destructor
    virtual ~YoungModulusHofmanD9();


    // Member Functions
    
        //- Update conductivity  
        virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr) ;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
