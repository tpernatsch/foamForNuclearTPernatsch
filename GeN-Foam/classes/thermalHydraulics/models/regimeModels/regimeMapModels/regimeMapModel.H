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
    Foam::regimeMapModel

Description
    Run-time selectable class to handle flow regime map models. It consists 
    of a HashTable of regimes indexed by regime name and a HashTable of
    volScalarFields indexed by pair name. These fields are the global
    field representing the drag coefficients (Kd), heat transfer coefficient
    (htc), interfacial area (iA), and are indexed in the table by the pair
    name. Each of the fields is obtained by combining the corresponding fields
    computed by the corresponding models within each regime object, depending
    on the regime distribution.
    The functionality of this class lies in:
    -   correctMap(), a pure abstract to be implemented in a spetialization
        of this class. Its main purpose is to adjust the scalarField of each
        regime (each regime is a scalarField, see regime.H for more info) so to
        represent the current regime distribution within the domain, for each
        regime;
    -   correctFieldTables(), a concrete implementation that adjusts the 
        volScalarField tables depending on the regime "distrbution".

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    regimeMapModel.C

\*---------------------------------------------------------------------------*/

#ifndef regimeMapModel_H
#define regimeMapModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "HashPtrTable.H"
#include "runTimeSelectionTables.H"
#include "stringOps.H"
#include "regime.H"
#include "typedefs.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class regimeMapModel Declaration
\*---------------------------------------------------------------------------*/

class regimeMapModel
    : public IOdictionary
{
protected :

    //- Ref to mesh
    const fvMesh& mesh_;

    //- Ref to dictionary in which all the physics models are specified on a
    //  regime-by-regime basis. Required as it is regimeMap that creates the
    //  regimes.
    const dictionary& physicsModelsDict_;

    //- Table of regimes indexed by regime name
    regimeTable regimes_;

public:

    //- Runtime type information
    TypeName("regimeMapModel");

    //- Runtime selection table declaration
    declareRunTimeSelectionTable
    (
        autoPtr,
        regimeMapModel,
        regimeMapModels,
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& physicsModelsDict
        ),
        (mesh, dict, physicsModelsDict)
    );

    //- Constructors

        regimeMapModel
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& physicsModelsDict
        );

        //- Destructor
        virtual ~regimeMapModel(){}

    //- Selectors
    
        //- The selector does not have to be the same as the constructor
        static autoPtr<regimeMapModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& physicsModelsDict
        );

    //- Members

        //- Updated the marker field of each regime in regimeTable
        //  regimeTable according to the desired map model
        virtual void correct() = 0;

        //- Needs to be called at the end of every correct. If a
        //  non-interpolated regime is present, set its requiresModelCorrection
        //  flag to true. If an interpolated regime is present, set the
        //  the requiresModelCorrection flag of the regimes from which it is
        //  interpolated from to true
        void setRequiresModelCorrection();

        //- Access
            regimeTable& regimes()
            {
                return regimes_;
            }

            const regimeTable& regimes() const
            {
                return regimes_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
