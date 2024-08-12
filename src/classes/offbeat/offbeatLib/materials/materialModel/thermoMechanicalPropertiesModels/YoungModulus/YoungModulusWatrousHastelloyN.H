
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
    Foam::YoungModulusWatrousHastelloyN

Description
    Class modelling Young Modulus of HN from work of Watrous.

SourceFiles
    YoungModulusWatrousHastelloyN.C

 \mainauthor
    Alejandra de Lara - University of Cambridge

 \contribution
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour)\n

 \date
    Sep 2023


\*---------------------------------------------------------------------------*/

#ifndef YoungModulusWatrousHastelloyN_H
#define YoungModulusWatrousHastelloyN_H

#include "YoungModulusModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class YoungModulusWatrousHastelloyN Declaration
\*---------------------------------------------------------------------------*/

class YoungModulusWatrousHastelloyN
:
    public YoungModulusModel
{
    // Private data



        //- Parameters for the Matpro model
        scalar par1;
        scalar par2;
        scalar par3;
        scalar par4;
   
        //- Perturbation parameter
        scalar perturb;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        YoungModulusWatrousHastelloyN(const YoungModulusWatrousHastelloyN&);

        //- Disallow default bitwise assignment
        void operator=(const YoungModulusWatrousHastelloyN&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("WatrousHastelloyN");

    // Constructors

        //- Construct from dictionary
        YoungModulusWatrousHastelloyN
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word defaultModel
        );
    //- Destructor
    virtual ~YoungModulusWatrousHastelloyN();


    // Member Functions
    
    //- Update conductivity
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr) ;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
