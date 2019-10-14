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
    Foam::dragModel

Description
    Class to handle interfacial drag between either the two fluid phases
    or a fluid phase and the structure. It is run-time selectable depending
    on whether a fluid-fluid (FF in short) or a fluid-structure (FS in short)
    pair is provided. For this reason, the class has two constructors
    and, consequently, two selector implementations in newDragModel.C. However,
    as the return type of Kd is always a tensorField, only one base class was
    deemed necessary, in order to keep it simple.
    Fluid-fluid models are to be added to the FF selection table, while 
    fluid-structure models to the FS selection table.

    By Stefan Radman (stefan.radman@epfl.ch)

SourceFiles
    dragModel.C
    newDragModel.C

\*---------------------------------------------------------------------------*/

#ifndef dragModel_H
#define dragModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluid.H"
#include "structureModel.H"
#include "FFPair.H"
#include "FSPair.H"
#include "zeroGradientFvPatchFields.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dragModel Declaration
\*---------------------------------------------------------------------------*/

class dragModel
:
    public IOdictionary
{
protected:

    //- Ref to mesh
    const fvMesh& mesh_;

    //- Name of phase pair
    const word pairName_;

    //- Name of first phase in pair
    const word name1_;

    //- Name of second phase in pair
    const word name2_;

    //- True if the drag models is between a fluid and structure region
    const bool withStructure_;

public:

    TypeName("dragModel");

    //- Static data members
    static const dimensionSet dimK;

    //- Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            dragModel,
            FFDragModels,
            (
                const objectRegistry& objReg,
                const dictionary& dragModelsDict,
                const FFPair& FFPair
            ),
            (objReg, dragModelsDict, FFPair)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            dragModel,
            FSDragModels,
            (
                const objectRegistry& objReg,
                const dictionary& dragModelsDict,
                const FSPair& FSPair,
                const word region
            ),
            (objReg, dragModelsDict, FSPair, region)
        );

    //- Constructors

        dragModel
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FFPair& FFPair
        );

        dragModel
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FSPair& FSPair,
            const word region
        );

        //- Destructor
        virtual ~dragModel(){}

    //- Selectors

        static autoPtr<dragModel> New
        (
            const objectRegistry& objReg,
            const dictionary& dragModelsDict,
            const FFPair& FFPair
        );

        static autoPtr<dragModel> New
        (
            const objectRegistry& objReg,
            const dictionary& dragModelsDict,
            const FSPair& FSPair,
            const word region
        );

        //- Make in table given a key, dict, source and destination object
        //  registries. Source registers are those where the FF/FS pairs will
        //  be looked for. Destination ones are those in which the model 
        //  IOobject will be registered
        static void makeInTable
        (
            const objectRegistry& srcObjReg,
            const objectRegistry& dstObjReg,
            const word key,
            const dictionary& dict,
            HashTable<autoPtr<dragModel>, word, word::hash>& table
        );

    //- Member Functions

        //- Drag coefficient used in momentum equations
        //    ddt(alpha1*rho1*U1) + ... = ... Kd&(U1-U2)
        //    ddt(alpha2*rho2*U2) + ... = ... Kd&(U2-U1)
        virtual tmp<volTensorField> Kd() const = 0;

        //- Access

            const word pairName() const
            {
                return pairName_;
            }

            const word name1() const
            {
                return name1_;
            }

            const word name2() const
            {
                return name2_;
            }

            bool withStructure() const
            {
                return withStructure_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
