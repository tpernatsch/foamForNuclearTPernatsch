/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::thermalExpansionConstant

Description
    Model for constant thermal expansion coefficient. The value is read from 
    dictionary.

SourceFiles
    thermalExpansionConstant.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef thermalExpansionConstant_H
#define thermalExpansionConstant_H

#include "thermalExpansionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermalExpansionConstant Declaration
\*---------------------------------------------------------------------------*/

class thermalExpansionConstant
:
    public thermalExpansionModel
{
    // Private data 
    
        //- Constant linear thermal expansion coefficient [1/K]
        scalar alpha_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermalExpansionConstant(const thermalExpansionConstant&);

        //- Disallow default bitwise assignment
        void operator=(const thermalExpansionConstant&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("constant");

    // Constructors

        //- Construct from dictionary
        thermalExpansionConstant
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict, 
            const word defaultModel   
        );

    //- Destructor
    virtual ~thermalExpansionConstant();


    // Member Functions
    
    //- Update thermalExpansion  
    virtual void correct
    (symmTensorField& sf, const scalarField& T, const labelList& addr) const;   
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
