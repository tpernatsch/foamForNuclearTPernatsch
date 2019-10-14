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
    Foam::heatTransferModel

Description
    Class to handle interfacial heat transfer between either the two fluid 
    phases or a fluid phase and the structure. It is run-time selectable 
    depending on whether a fluid-fluid (FF in short) or a fluid-structure 
    (FS in short) model is to be selected. For this reason, the class has two 
    constructors and, consequently, two selector implementations in 
    newHeatTransferModel.C
    Fluid-fluid models are to be added to the FF selection table, while 
    fluid-structure models to the FS selection table.

    By Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    heatTransferModel.C
    newHeatTransferModel.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferModel_H
#define heatTransferModel_H

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
                           Class heatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class heatTransferModel
:
    public IOdictionary
{
protected:

        //- Ref to mesh
    const fvMesh& mesh_;

    //- Name of phase pair
    const word pairName_;

    //- Name of the bulk phase
    const word nameBulk_;

    //- Name of the phase across the interface
    const word nameInterface_;

    //- True if the drag models is between a fluid and structure region
    const bool withStructure_;

public:

    TypeName("heatTransferModel");

    //- Static data members
    static const dimensionSet dimHtc;

    //- Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            heatTransferModel,
            FFHeatTransferModels,
            (
                const objectRegistry& objReg,
                const dictionary& heatTransferModelsDict,
                const FFPair& FFPair,
                const word nameBulk,
                const word nameInterface
            ),
            (objReg, heatTransferModelsDict, FFPair, nameBulk, nameInterface)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            heatTransferModel,
            FSHeatTransferModels,
            (
                const objectRegistry& objReg,
                const dictionary& heatTransferModelsDict,
                const FSPair& FSPair,
                const word region
            ),
            (objReg, heatTransferModelsDict, FSPair, region)
        );

    //- Constructors

        heatTransferModel
        (
            const objectRegistry& objReg,
            const dictionary& heatTransferModelsDict,
            const FFPair& FFPair,
            const word nameBulk,
            const word nameInterface
        );

        heatTransferModel
        (
            const objectRegistry& objReg,
            const dictionary& heatTransferModelsDict,
            const FSPair& FSPair,
            const word region
        );

        //- Destructor
        virtual ~heatTransferModel(){}

    //- Selectors

        static autoPtr<heatTransferModel> New
        (
            const objectRegistry& objReg,
            const dictionary& heatTransferModelsDict,
            const FFPair& FFPair,
            const word nameBulk,
            const word nameInterface
        );

        static autoPtr<heatTransferModel> New
        (
            const objectRegistry& objReg,
            const dictionary& heatTransferModelsDict,
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
            HashTable<autoPtr<heatTransferModel>, word, word::hash>& table
        );

    //- Member Functions

        //- Heat transfer coefficient used in the energy equations
        //  (iA12() is the interfacial area)
        //    ddt(alpha1*rho1*h1) + ... = ... iA12()*htc12()*(T1-T2)
        virtual tmp<volScalarField> htc() const = 0;

    //- Access

        const word pairName() const
        {
            return pairName_;
        }

        const word nameBulk() const
        {
            return nameBulk_;
        }

        const word nameInterface() const
        {
            return nameInterface_;
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
