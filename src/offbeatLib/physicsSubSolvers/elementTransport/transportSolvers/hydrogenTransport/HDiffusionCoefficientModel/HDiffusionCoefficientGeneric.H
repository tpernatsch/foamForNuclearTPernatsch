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
    Foam::creepModel

Description


Usage

\mainauthor
    Matthieu Reymond - EPFL

\contribution

\date 
    Feb. 2023

\*---------------------------------------------------------------------------*/

#ifndef HDiffusionCoefficientGeneric_H
#define HDiffusionCoefficientGeneric_H

#include "HDiffusionCoefficientModel.H"
#include "fvMesh.H"
#include "Table.H"
#include "zeroGradientFvPatchField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class "HDiffusionCoefficientGeneric" Declaration
\*---------------------------------------------------------------------------*/

class HDiffusionCoefficientGeneric
:
    public HDiffusionCoefficientModel
{

protected:
        
        scalar D0_;
        scalar Ed_;
        
        //- Disallow default bitwise copy construct
        HDiffusionCoefficientGeneric(const HDiffusionCoefficientGeneric&);

        //- Disallow default bitwise assignment
        void operator=(const HDiffusionCoefficientGeneric&);
        

public:
 
    //- Runtime type information
    TypeName("HDiffusionCoefficientGeneric");

    // Constructors

        //- Construct from dictionary
        HDiffusionCoefficientGeneric
        (
            const fvMesh& mesh,
            const materials& mat,
            const dictionary& lawDict
        );

    //- Destructor
   virtual ~HDiffusionCoefficientGeneric();

    // Member Functions    
    
    //- Update TSSd
    virtual void correctDH
    (
        const labelList& addr
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
