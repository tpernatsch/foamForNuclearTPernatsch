/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::byRegime

Description

Author
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef byRegimeModel_H
#define byRegimeModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regimeMapModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class byRegime Declaration
\*---------------------------------------------------------------------------*/

template<class valueType, class baseModel>
class byRegimeModel
:
    public baseModel
{
protected:

    //- Ref to regime map
    const regimeMapModel& regimeMap_;

    //- List of twoPhaseDragMultiplier models indexed by regime label
    List<autoPtr<baseModel>> modelPtrs_;

    //- Local object registry in which sub-models (i.e. models in modelPtrs_
    //  are registered to avoid name conflicts if all registered in the mesh
    //  object registry (yeah, but do I reaaally need this? Eh!))
    objectRegistry modelObjReg_;

public:

    //- Runtime type information
    //TypeName("byRegime");

    // Constructors

        byRegimeModel
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        template<class firstArgType>
        byRegimeModel
        (
            const firstArgType& firstArg,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~byRegimeModel(){}

    // Member Functions

        //- Return value predicted by model
        virtual valueType value(const label& celli) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "byRegimeModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
