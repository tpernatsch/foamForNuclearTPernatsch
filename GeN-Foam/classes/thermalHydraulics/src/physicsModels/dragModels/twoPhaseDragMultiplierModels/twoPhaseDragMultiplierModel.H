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
    Foam::twoPhaseDragMultiplierModel

Description

Author
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    twoPhaseDragMultiplierModel.C
    newTwoPhaseDragMultiplierModel.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseDragMultiplierModel_H
#define twoPhaseDragMultiplierModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "fluid.H"
#include "FSPair.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class twoPhaseDragMultiplierModel Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseDragMultiplierModel
:
    public IOdictionary
{
protected:

    //- Ref to mesh
    const fvMesh& mesh_;

    //- Ptr to fluid whose properties are used as basis for the calculation of
    //  the two phase pressure drop
    const fluid* mFluidPtr_;

    //- Ptr to other fluid in the system
    const fluid* oFluidPtr_;

    //- Ptr to multiplier fluid (single-phase) friction factor tensor
    const volTensorField* mKd_;

    //- Maximum allowable value of phi2
    scalar maxPhi2_;

public:

    //- Runtime type information
    TypeName("twoPhaseDragMultiplierModel");

    //- Declare runtime construction
        
        declareRunTimeSelectionTable
        (
            autoPtr,
            twoPhaseDragMultiplierModel,
            twoPhaseDragMultiplierModels,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                const objectRegistry& objReg
            ),
            (mesh, dict, objReg)
        );

    //- Constructors

        twoPhaseDragMultiplierModel
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~twoPhaseDragMultiplierModel(){}

    //- Selectors

        static autoPtr<twoPhaseDragMultiplierModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const objectRegistry& objReg
        );

    //- Member Functions

        //- True if flow is essentially single phase in celli
        bool onePhase(const label& celli) const;

        //- Return multiplier
        virtual scalar phi2(const label& celli) const;

        //- Return total drag factor tensor divided by
        //  multiplier fluid velocity.
        virtual tensor value
        (
            const label& celli
        ) const;

        //-
        void correctField(volTensorField& KdTotByU);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
