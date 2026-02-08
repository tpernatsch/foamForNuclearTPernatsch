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
    Foam::RoutbortFastMOXCreepModel 

Description
    Creep model for fast-MOX fuel.

    Routbort ref : https://doi.org/10.1016/0022-3115(72)90036-0
    Combette ref : https://doi.org/10.1016/0022-3115(77)90039-3

SourceFiles 
    RoutbortFastMOXCreepModel.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro - EPFL\n

\date 
    July 2024

\*---------------------------------------------------------------------------*/

#ifndef RoutbortFastMOXCreepModel_H
#define RoutbortFastMOXCreepModel_H

#include "creepModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam 
{

/*---------------------------------------------------------------------------*\
                    Class "RoutbortFastMOXCreepModel" Declaration
\*---------------------------------------------------------------------------*/

class RoutbortFastMOXCreepModel
:
    public creepModel
{
    // Private Data

        //- Irradiation creep Strain 
        volSymmTensorField& epsilonCreepIrr_;

        //- Thermal creep Strain 
        volSymmTensorField& epsilonCreepTh_;
        
        //- equivalent Irradiation Creep Strain 
        volScalarField& epsilonCreepIrrEq_;
        
        //- Increment of equivalent Irradiation Creep Strain 
        volScalarField& DepsilonCreepIrrEq_;
        
        //- Equivalent Secondary Thermal Creep Strain 
        volScalarField& epsilonCreepThEq_;
        
        //- Increment of equivalent Secondary Thermal Creep Strain 
        volScalarField& DepsilonCreepThEq_;
        
        //- Name heatSource field
        const word heatSourceName_;
        
        //- Name grain size field
        const word grainRadiusName_;        

        //- Model variables
        const scalar A_;
        const scalar B_;
        const scalar C_;
        const scalar Q1_;
        const scalar Q2_;
        const scalar R_;

        //- Energy of fission event
        const scalar Efission_;

        //- Switch for irradiation creep
        Switch irradiationCreep_;

        //- Index of current time
        scalar timeIndex_;

        //- Relaxation factor
        scalar relax_;

        //- Activate Newton-Raphson method for convergence
        Switch NewtonRaphsonMethod_;

        //- Tolerance for inner convergence of the Newton-Raphson method
        scalar NewtonRaphsonTolerance_;

        const bool creepRateLimiter_;
        const scalar creepRateLimit_;
        
    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        RoutbortFastMOXCreepModel(const RoutbortFastMOXCreepModel&);

        //- Disallow default bitwise assignment
        void operator=(const RoutbortFastMOXCreepModel&);

public:
 
    //- Runtime type information
    TypeName("UPuO2Routbort");


    // Constructors

        //- Construct from mesh and thermo-mechanical properties
        RoutbortFastMOXCreepModel
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );


    //- Destructor
    ~RoutbortFastMOXCreepModel();


    // Member Functions

        //- Update mechanical properties according to rheology model        
        virtual void correctCreep
        (
            volSymmTensorField& epsilonEl,
            const labelList& addr
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //