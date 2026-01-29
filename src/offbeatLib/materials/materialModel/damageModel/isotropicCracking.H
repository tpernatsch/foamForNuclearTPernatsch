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
    Foam::isotropicCracking

Description
Usage
    \verbatim
    rheology byMaterial;
    
    ...

    // List of materials, one per cellZone.
    materials
    {
        ...

        ...
    }
    \endverbatim

SourceFiles
    isotropicCracking.c

\todo 

\mainauthor
    Alessandro Scolaro (LRS -- EPFL)

\contribution

\date 
    March 2024.

\*---------------------------------------------------------------------------*/
#ifndef EXPONENTIALSOFTENINGCRACKING_H
#define EXPONENTIALSOFTENINGCRACKING_H

#include "damageModel.H"
#include "Table.H"
#include "Time.H"
#include "mechanicsSubSolver.H"
#include "sliceMapper.H"

namespace Foam{

/*---------------------------------------------------------------------------*\
                         Class isotropicCracking Declaration
\*---------------------------------------------------------------------------*/


class isotropicCracking
:
    public damageModel
{
    // Private Member Functions        

    //- Disallow default bitwise copy construct
    isotropicCracking(const isotropicCracking&);

    //- Disallow default bitwise assignment
    void operator=(const isotropicCracking&);

protected:

    //- Pointer to mapper class
    const sliceMapper* mapper_;

    //- Reference to nCracks field
    volScalarField& nCracks_;
        
    //- Decide between alternative isotropic cracking models 
    word isotropicCrackingType_;

    //- Maximum number of cracks (from 0 to 12) for Barani model
    scalar nCracksMax_;

    //- Fixed number of cracks for JankusWeeks model
    scalar nCracksFixed_;

public:
    //- Runtime type information
    TypeName("isotropicCracking");

    // Constructors

        //- Construct from mesh, global options, materials , dict and labelList
        isotropicCracking
        (
            const fvMesh& mesh,
            const dictionary& damageDict,
            const dictionary& baseDict
        );

    //- Destructor
    ~isotropicCracking();


    // Member Functions

        //- Soften E and nu
        virtual void soften
        (
            scalarField& E,
            scalarField& nu,
            const labelList& addr
        );
};

} // End namespace Foam

#endif
