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
    Foam::heatCapacityMolybdenum

Description

SourceFiles
    heatCapacityMolybdenum.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef heatCapacityMolybdenum_H
#define heatCapacityMolybdenum_H

#include "heatCapacityModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heatCapacityMolybdenum Declaration
\*---------------------------------------------------------------------------*/

class heatCapacityMolybdenum
:
    public heatCapacityModel
{
    // Private data

        //- Parameters for the model
        scalar par1;
        scalar par2;
        scalar par3;

        //- Perturbation parameter
        scalar perturb;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        heatCapacityMolybdenum(const heatCapacityMolybdenum&);

        //- Disallow default bitwise assignment
        void operator=(const heatCapacityMolybdenum&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("Mo");

    // Declare run-time constructor selection table

    // Constructors

        //- Construct from dictionary
        heatCapacityMolybdenum
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel   
        );


    // Selectors

    //- Destructor
    virtual ~heatCapacityMolybdenum();


    // Member Functions
    
    //- Update heatCapacity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
