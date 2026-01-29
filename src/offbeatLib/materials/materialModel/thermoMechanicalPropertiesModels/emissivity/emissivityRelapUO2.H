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
    Foam::emissivityRelapUO2

Description
    UO2 emissivity model derived from Relap.

SourceFiles
    emissivityRelapUO2.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef emissivityRelapUO2_H
#define emissivityRelapUO2_H

#include "emissivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class emissivityRelapUO2 Declaration
\*---------------------------------------------------------------------------*/

class emissivityRelapUO2
:
    public emissivityModel
{
    // Private data

        //- Parameters for the Relap model
        scalar par1;
        scalar par2;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        emissivityRelapUO2(const emissivityRelapUO2&);

        //- Disallow default bitwise assignment
        void operator=(const emissivityRelapUO2&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UO2Relap");

    // Constructors

        //- Construct from dictionary
        emissivityRelapUO2
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict, 
            const word defaultModel   
        );
        
    //- Destructor
    virtual ~emissivityRelapUO2();


    // Member Functions
    
    //- Update emissivity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
