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
    Foam::yieldStressConstant

Description
    User defined constant yield stress.

    \htmlinclude yieldStressModel_yieldStressConstant.html

SourceFiles
    yieldStressConstant.C

\mainauthor
    Matthieu Reymond

\contribution
    
\date
    Feb. 2023 


\*---------------------------------------------------------------------------*/

#ifndef yieldStressconstant_H
#define yieldStressConstant_H

#include "yieldStressModel.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class "yieldStressConstant" Declaration
\*---------------------------------------------------------------------------*/

class yieldStressConstant
:
    public yieldStressModel 
{    

    // Private data    
    // Private Member Functions
        
        //- Value of constant yieldStress [Pa]
        scalar sigmaYValue_;

        //- Disallow default bitwise copy construct
        yieldStressConstant(const yieldStressConstant&);

        //- Disallow default bitwise assignment
        void operator=(const yieldStressConstant&);
        
protected:
        
        // Protected data
    
public:
 
    //- Runtime type information
    TypeName("yieldStressConstant");

    // Constructors

        //- Construct from dictionary
        yieldStressConstant
        (
            const fvMesh& mesh,
            const dictionary& lawDdict
        );

    //- Destructor
   virtual ~yieldStressConstant();

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
