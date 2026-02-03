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
    Foam::conductivityUO2IFA601

Description
    Model for conductivity of UO2 - IFA-601 rod.

SourceFiles
    conductivityUO2IFA601.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef conductivityUO2IFA601_H
#define conductivityUO2IFA601_H

#include "conductivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class conductivityUO2IFA601 Declaration
\*---------------------------------------------------------------------------*/

class conductivityUO2IFA601
:
    public conductivityModel
{
    // Private data
        
        //- Name used for the burnup field
        const word burnupName_;
        
        //- Reference to Burnup field
        const volScalarField* Bu_;

        //- Parameters for the Matpro model
        scalar par1;
        scalar par2;
        scalar par3;
        scalar par4;
        scalar par5;
        scalar par6;
        scalar par7;
        scalar par8;
        
    // Private Member Functions

        //- Disallow default bitwise copy construct
        conductivityUO2IFA601(const conductivityUO2IFA601&);

        //- Disallow default bitwise assignment
        void operator=(const conductivityUO2IFA601&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UO2Ifa601");

    // Constructors

        //- Construct from dictionary
        conductivityUO2IFA601
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict,
            const word defaultModel         
        );
        
    //- Destructor
    virtual ~conductivityUO2IFA601();


    // Member Functions
    
    //- Update conductivity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
