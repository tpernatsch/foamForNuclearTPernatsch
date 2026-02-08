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
    Foam::heatCapacityMatproUO2

Description
    Heat Capacity model for UO2 derived from Matprov11.

SourceFiles
    heatCapacityMatproUO2.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef heatCapacityMatproUO2_H
#define heatCapacityMatproUO2_H

#include "heatCapacityModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heatCapacityMatproUO2 Declaration
\*---------------------------------------------------------------------------*/

class heatCapacityMatproUO2
:
    public heatCapacityModel
{
    // Private data

        //- Parameters for the Matpro model
        scalar OM;
        scalar K1;
        scalar K2;
        scalar K3;
        scalar phi;
        scalar ED;
        scalar R;
        scalar par1;
        scalar par2;
        scalar par3;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        heatCapacityMatproUO2(const heatCapacityMatproUO2&);

        //- Disallow default bitwise assignment
        void operator=(const heatCapacityMatproUO2&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UO2Matpro");

    // Declare run-time constructor selection table

    // Constructors

        //- Construct from dictionary
        heatCapacityMatproUO2
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict, 
            const word defaultModel   
        );


    // Selectors

    //- Destructor
    virtual ~heatCapacityMatproUO2();


    // Member Functions
    
    //- Update heatCapacity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
