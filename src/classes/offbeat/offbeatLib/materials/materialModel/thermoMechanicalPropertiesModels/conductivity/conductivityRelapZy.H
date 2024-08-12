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
    Foam::conductivityRelapZy

Description
    Model for conductivity of Zircaloy material derived by RELAP.
    
SourceFiles
    conductivityRelapZy.C

\todo
    TODO: change name to MATPRO

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef conductivityRelapZy_H
#define conductivityRelapZy_H

#include "conductivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class conductivityRelapZy Declaration
\*---------------------------------------------------------------------------*/

class conductivityRelapZy
:
    public conductivityModel
{
    // Private data

        //- Parameters for the Relap model
        scalar par1;
        scalar par2;
        scalar par3;
        scalar par4;
        
        //- Perturbation factor
        scalar perturb;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        conductivityRelapZy(const conductivityRelapZy&);

        //- Disallow default bitwise assignment
        void operator=(const conductivityRelapZy&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("ZyRELAP");

    // Constructors

        //- Construct from dictionary
        conductivityRelapZy
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel   
        );

    //- Destructor
    virtual ~conductivityRelapZy();


    // Member Functions
    
    //- Update conductivity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
