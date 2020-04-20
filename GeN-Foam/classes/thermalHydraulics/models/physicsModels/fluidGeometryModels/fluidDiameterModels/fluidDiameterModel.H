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
    Foam::fluidDiameterModel

Description
    Class to handle the charachteristic dimension of a dispersed fluid phase. 
    This is used to compute the two-phase Reynolds number and fluid-fluid
    heat transfer coefficients if correlations are based on the Nusselt 
    number. It is alway limited by the hydraulic diameter of the structure the 
    fluid is in, if applicable.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    fluidDiameterModel.C
    newFluidDiameterModel.C

\*---------------------------------------------------------------------------*/

#ifndef fluidDiameterModel_H
#define fluidDiameterModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "zeroGradientFvPatchFields.H"
#include "fluid.H"
#include "FFPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fluidDiameterModel Declaration
\*---------------------------------------------------------------------------*/

class fluidDiameterModel
:
    public IOdictionary
{
protected:

    //- Ref to fluid
    const fvMesh& mesh_;

    //- Ref to fluid
    const fluid& fluid_;

public:

    //- Runtime type information
    TypeName("fluidDiameterModel");

    //- Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluidDiameterModel,
            fluidDiameterModels,
            (
                const objectRegistry& objReg,
                const dictionary& dict,
                const fluid& fluid
            ),
            (objReg, dict, fluid)
        );

    //- Constructors

        fluidDiameterModel
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fluid& fluid
        );

        //- Destructor
        virtual ~fluidDiameterModel(){}

    //- Selectors

        static autoPtr<fluidDiameterModel> New
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fluid& fluid
        );

    //- Member Functions

        //- Fluid diameter
        virtual tmp<volScalarField> Dh() const = 0;

        //- Access
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
