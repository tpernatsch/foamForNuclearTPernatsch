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
    Foam::emissivityConstant

Description
    Model for constant emissivity. The value is read from dictionary.

SourceFiles
    emissivityConstant.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef emissivityConstant_H
#define emissivityConstant_H

#include "emissivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class emissivityConstant Declaration
\*---------------------------------------------------------------------------*/

class emissivityConstant
:
    public emissivityModel
{
    // Private data
        
        //- Constant emmissivity [-]
        scalar emissivity_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        emissivityConstant(const emissivityConstant&);

        //- Disallow default bitwise assignment
        void operator=(const emissivityConstant&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("constant");

    // Constructors

        //- Construct from dictionary
        emissivityConstant
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict, 
            const word defaultModel   
        );
        
    //- Destructor
        virtual ~emissivityConstant();


    // Member Functions
    
    //- Update emissivity  
        virtual void correct
        (scalarField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
