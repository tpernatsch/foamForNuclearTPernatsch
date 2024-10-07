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
    Foam::heatCapacityFinkUPuO2

Description
    Heat Capacity model for (U,Pu)O2 derived from Fink.
    Source : https://info.ornl.gov/sites/publications/Files/Pub57523.pdf

SourceFiles
    heatCapacityFinkUPuO2.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef heatCapacityFinkUPuO2_H
#define heatCapacityFinkUPuO2_H

#include "heatCapacityModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heatCapacityFinkUPuO2 Declaration
\*---------------------------------------------------------------------------*/

class heatCapacityFinkUPuO2
:
    public heatCapacityModel
{
    // Private data

        //- Parameters for the Fink PuO2 heat capacity model
        scalar C1_;
        scalar C2_;
        scalar C3_;
        scalar theta_;
        scalar Ea_;

        //- Perturbation parameter
        scalar perturb;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        heatCapacityFinkUPuO2(const heatCapacityFinkUPuO2&);

        //- Disallow default bitwise assignment
        void operator=(const heatCapacityFinkUPuO2&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UPuO2Fink");

    // Declare run-time constructor selection table

    // Constructors

        //- Construct from dictionary
        heatCapacityFinkUPuO2
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel   
        );


    // Selectors

    //- Destructor
    virtual ~heatCapacityFinkUPuO2();


    // Member Functions
    
    //- Update heatCapacity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
