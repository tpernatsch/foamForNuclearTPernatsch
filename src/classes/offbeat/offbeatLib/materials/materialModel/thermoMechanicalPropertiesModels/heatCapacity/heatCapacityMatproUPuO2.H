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
    Foam::heatCapacityMatproUPuO2

Description
    Heat Capacity model for (U,Pu)O2 derived from MATPROv11.

SourceFiles
    heatCapacityMatproUPuO2.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef heatCapacityMatproUPuO2_H
#define heatCapacityMatproUPuO2_H

#include "heatCapacityModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heatCapacityMatproUPuO2 Declaration
\*---------------------------------------------------------------------------*/

class heatCapacityMatproUPuO2
:
    public heatCapacityModel
{
    // Private data

        //- Parameters for the Matpro PuO2 heat capacity model
        scalar OM_;
        scalar K1_;
        scalar K2_;
        scalar K3_;
        scalar theta_;
        scalar ED_;
        scalar R_;

        //- Weight fractions of UO2 and PuO2 in the MOX
        scalar wUO2_;
        scalar wPuO2_;


        //- Perturbation parameter
        scalar perturb;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        heatCapacityMatproUPuO2(const heatCapacityMatproUPuO2&);

        //- Disallow default bitwise assignment
        void operator=(const heatCapacityMatproUPuO2&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UPuO2MATPRO");

    // Declare run-time constructor selection table

    // Constructors

        //- Construct from dictionary
        heatCapacityMatproUPuO2
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel   
        );


    // Selectors

    //- Destructor
    virtual ~heatCapacityMatproUPuO2();


    // Member Functions
    
    //- Update heatCapacity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
