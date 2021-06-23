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
    Foam::dispersionModel

Description
    Class to handle fluid-structure interfacial contact contact fractions

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    dispersionModel.C
    newdispersionModel.C

\*---------------------------------------------------------------------------*/

#ifndef dispersionModel_H
#define dispersionModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "fluid.H"

namespace Foam
{

class FFPair;

/*---------------------------------------------------------------------------*\
            Class dispersionModel Declaration
\*---------------------------------------------------------------------------*/

class dispersionModel
:
    public IOdictionary
{
protected:

    const fvMesh& mesh_;

    const FFPair& pair_;

public:

    //- Runtime type information
    TypeName("dispersionModel");

    //- Declare runtime construction
        
        declareRunTimeSelectionTable
        (
            autoPtr,
            dispersionModel,
            dispersionModels,
            (
                const FFPair& pair,
                const dictionary& dict,
                const objectRegistry& objReg
            ),
            (pair, dict, objReg)
        );

    // Constructors

        dispersionModel
        (
            const FFPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~dispersionModel();

    // Selectors

        static autoPtr<dispersionModel> New
        (
            const FFPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

    // Member Functions

        //- Return dispersion value of fluid1_ (ALWAYS fluid1)
        virtual scalar value(const label& celli) const = 0;

        void correctField(scalarField& d);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
