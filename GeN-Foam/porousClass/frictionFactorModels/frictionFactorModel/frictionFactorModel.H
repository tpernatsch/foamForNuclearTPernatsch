/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::frictionFactorModel

Description

SourceFiles
    frictionFactorModel.C
    newfrictionFactorModel.C

\*---------------------------------------------------------------------------*/

#ifndef frictionFactorModel_H
#define frictionFactorModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regIOobject.H"
#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class frictionFactorModel Declaration
\*---------------------------------------------------------------------------*/

class frictionFactorModel
:
    public regIOobject
{

public:

    //- Runtime type information
    TypeName("frictionFactorModel");

    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        frictionFactorModel,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (dict, mesh)
    );

    // Static data members
    static const dimensionSet dimK;

    // Constructors
    frictionFactorModel
    (
        const dictionary& dict,
        const fvMesh& mesh
    );

    //- Destructor
    virtual ~frictionFactorModel();

    // Selectors
    static autoPtr<frictionFactorModel> New
    (
        const dictionary& dict,
        const fvMesh& mesh
    );

    // Member Functions
    virtual tmp<volScalarField> frictionFactor() const = 0;
    bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
