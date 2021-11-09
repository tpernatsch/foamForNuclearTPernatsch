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
    Foam::virtualMassCoefficientModel

Description
    Class to handle fluid-fluid interfacial area

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    virtualMassCoefficientModel.C
    newInterfacialAreaModel.C

\*---------------------------------------------------------------------------*/

#ifndef virtualMassCoefficientModel_H
#define virtualMassCoefficientModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluid.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class FFPair;

/*---------------------------------------------------------------------------*\
            Class virtualMassCoefficientModel Declaration
\*---------------------------------------------------------------------------*/

class virtualMassCoefficientModel
:
    public IOdictionary
{
protected:

    const fvMesh& mesh_;

    const FFPair& pair_;

    const fluid& fluid1_;

    const fluid& fluid2_;

public:

    //- Runtime type information
    TypeName("virtualMassCoefficientModel");

    //- Declare runtime construction
        
        declareRunTimeSelectionTable
        (
            autoPtr,
            virtualMassCoefficientModel,
            virtualMassCoefficientModels,
            (
                const FFPair& pair,
                const dictionary& dict,
                const objectRegistry& objReg
            ),
            (pair, dict, objReg)
        );

    // Constructors

        virtualMassCoefficientModel
        (
            const FFPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~virtualMassCoefficientModel(){}

    // Selectors

        static autoPtr<virtualMassCoefficientModel> New
        (
            const FFPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

    // Member Functions

        //- Return fraction of structure interfacial area in
        //  contact with the fluid
        virtual scalar value(const label& celli) const = 0;

        //-
        virtual void correctField(volScalarField& Vm) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
