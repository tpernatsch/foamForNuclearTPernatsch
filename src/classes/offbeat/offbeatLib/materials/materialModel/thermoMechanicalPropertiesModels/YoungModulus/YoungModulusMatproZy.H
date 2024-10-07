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
    Foam::YoungModulusMatproZy

Description
    Class modelling Young Modulus of Zy from MATPROv11.

SourceFiles
    YoungModulusMatproZy.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef YoungModulusMatproZy_H
#define YoungModulusMatproZy_H

#include "YoungModulusModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class YoungModulusMatproZy Declaration
\*---------------------------------------------------------------------------*/

class YoungModulusMatproZy
:
    public YoungModulusModel
{
    // Private data

        //- Name of fast fluence field
        const word fastFluenceName_;

        //- Fast fluence field
        const volScalarField* phi_;

        //- Cold work
        scalar cW_;

        //- Average oxygen concentration minus as-received concentration (kg O2 / kg zircaloy) [-]
        scalar oxC_;

        //- Parameters for the Matpro model
        scalar par1;
        scalar par2;
        scalar par3;
        scalar par4;
        scalar par5;
        scalar par6;
        scalar par7;
        scalar par8;
        scalar par9;
        scalar par10;

        //- Perturbation parameter
        scalar perturb;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        YoungModulusMatproZy(const YoungModulusMatproZy&);

        //- Disallow default bitwise assignment
        void operator=(const YoungModulusMatproZy&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("ZyMATPRO");

    // Constructors

        //- Construct from dictionary
        YoungModulusMatproZy
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel   
        );
    //- Destructor
    virtual ~YoungModulusMatproZy();


    // Member Functions
    
    //- Update conductivity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr) ;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
