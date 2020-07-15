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
    Class to handle the charachteristic dimension of a dispersed fluid phase. 
    This is used to compute the two-phase Reynolds number and fluid-fluid
    heat transfer coefficients if correlations are based on the Nusselt 
    number. It is alway limited by the hydraulic diameter of the structure the 
    fluid is in, if applicable.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    twoPhaseDragMultiplierModel.C
    newtwoPhaseDragMultiplierModel.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseDragMultiplierModel_H
#define twoPhaseDragMultiplierModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "zeroGradientFvPatchFields.H"
#include "fluid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declaration
class fluidGeometry;

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

    //- Ref to multiplier fluid
    const fluid& mFluid_;

    //- Ptr to other fluid
    const fluid* oFluidPtr_;

    //- Two phase multiplier
    volScalarField phi2_;

    //- Maximum allowable value of phi2
    dimensionedScalar maxPhi2_;

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
                const objectRegistry& objReg,
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (objReg, dict, mesh)
        );

    //- Constructors

        twoPhaseDragMultiplierModel
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Destructor
        virtual ~twoPhaseDragMultiplierModel(){}

    //- Selectors

        static autoPtr<twoPhaseDragMultiplierModel> New
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Member Functions

        //- Limit phi2 by maxPhi
        void limitPhi2();

        //- Correct two phase multiplier phi2_
        virtual void correct() = 0;

        //- Correct Kd for a given fluid via two phase multiplier approach
        virtual void correctKdTable
        (
            const fluidGeometry& fG,
            volTensorFieldPtrTable& Kd
        ) const;

        //- Access
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
