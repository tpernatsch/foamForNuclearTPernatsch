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
    Foam::compoundFrictionFactorModel

Description

SourceFiles
    compoundFrictionFactorModel.C
    newcompoundFrictionFactorModel.C

\*---------------------------------------------------------------------------*/

#ifndef compoundFrictionFactorModel_H
#define compoundFrictionFactorModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#include "volFields.H"
#include "dictionary.H"
#include "frictionFactorModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class compoundFrictionFactorModel Declaration
\*---------------------------------------------------------------------------*/

class compoundFrictionFactorModel
:
    public regIOobject
{

protected:

    typedef PtrList<frictionFactorModel> frictionFactorModelList;

private:

    frictionFactorModelList frictionFactors_;

public:

    //- Runtime type information
    TypeName("compoundFrictionFactorModel");

    // Static data members
    static const dimensionSet dimK;

    // Constructor
    compoundFrictionFactorModel
    (
        const dictionary& dict,
        const fvMesh& mesh
    );

    //- Destructor
    ~compoundFrictionFactorModel();

    // Member Functions
    bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
