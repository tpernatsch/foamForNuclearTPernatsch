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
    Foam::hardening

Description
    User defined constant yield stress as a function of the equivalent plastic strain.

    \htmlinclude yieldStressModel_hardening.html

SourceFiles
    hardening.C

\mainauthor
    Matthieu Reymond - EPFL

\contribution

\date 
    Feb. 2023

\*---------------------------------------------------------------------------*/

#ifndef hardening_H
#define hardening_H

#include "yieldStressModel.H"
#include "fvMesh.H"
#include "Table.H"
#include "zeroGradientFvPatchField.H"

#ifdef OPENFOAMESI
    #include "Function1.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class "hardening" Declaration
\*---------------------------------------------------------------------------*/

class hardening
:
    public yieldStressModel
{
        //- Disallow default bitwise copy construct
        hardening(const hardening&);

        //- Disallow default bitwise assignment
        void operator=(const hardening&);
        
protected:
        
        // Protected data
#ifdef OPENFOAMFOUNDATION
        //- Table of post-yield stress versus plastic strain
        autoPtr<Function1s::Table<scalar>> stressPlasticStrainSeries_;
#elif OPENFOAMESI
        typedef Foam::Function1Types::Table<scalar> Table;

        //- Table of post-yield stress versus plastic strain
        autoPtr<Table> stressPlasticStrainSeries_;
#endif        
    
public:
 
    //- Runtime type information
    TypeName("hardening");

    // Constructors

        //- Construct from dictionary
        hardening
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );

    //- Destructor
   virtual ~hardening();

    // Member Functions    
    
        //- Update yield stress
        virtual void correctYieldStress
        (
            const labelList& addr
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
