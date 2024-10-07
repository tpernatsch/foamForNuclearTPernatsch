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
    Foam::heatCapacityMatproZy

Description
    Correlation for Zircaloy heat capacity from Matpro.
    Link : https://www.nrc.gov/docs/ML1429/ML14296A063.pdf - page 60.

SourceFiles
    heatCapacityMatproZy.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef heatCapacityMatproZy_H
#define heatCapacityMatproZy_H

#include "heatCapacityModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class heatCapacityMatproZy Declaration
\*---------------------------------------------------------------------------*/

class heatCapacityMatproZy
:
    public heatCapacityModel
{
    // Private data
    
        //- Perturbation parameter
        scalar perturb;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        heatCapacityMatproZy(const heatCapacityMatproZy&);

        //- Disallow default bitwise assignment
        void operator=(const heatCapacityMatproZy&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("ZyMATPRO");

    // Declare run-time constructor selection table

    // Constructors

        //- Construct from mesh, dicionary and defaultModel
        heatCapacityMatproZy
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel
        );


    // Selectors


    //- Destructor
    virtual ~heatCapacityMatproZy();


    // Member Functions
    
    //- Update heatCapacity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
