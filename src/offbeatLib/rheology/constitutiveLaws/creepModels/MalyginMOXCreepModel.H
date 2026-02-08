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
    Foam::MalyginMOXCreepModel 

Description
    Creep model for MOX fuel.

    Sources : 
    - https://link.springer.com/article/10.1007/s10512-010-9250-y
    - https://link.springer.com/article/10.1007/s10512-010-9266-3

SourceFiles 
    MalyginMOXCreepModel.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro - EPFL\n

\date 
    July 2024

\*---------------------------------------------------------------------------*/

#ifndef MalyginMOXCreepModel_H
#define MalyginMOXCreepModel_H

#include "creepModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam 
{

/*---------------------------------------------------------------------------*\
                    Class "MalyginMOXCreepModel" Declaration
\*---------------------------------------------------------------------------*/

class MalyginMOXCreepModel
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

        const word oxygenMetalRatioName_;
        const word youngModulusName_;
        const word porosityName_;     

        const scalar PuConc_;


        //- Model variables
        const scalar Efission_;
        const scalar KB_;
        const scalar DH_;
        const scalar Q_;
        const scalar Q1_;
        const scalar Q2_;
        const scalar Q3_;
        const scalar A1_;
        const scalar A2_;
        const scalar A3_;
        const scalar B1_;

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
        MalyginMOXCreepModel(const MalyginMOXCreepModel&);

        //- Disallow default bitwise assignment
        void operator=(const MalyginMOXCreepModel&);

public:
 
    //- Runtime type information
    TypeName("UPuO2Malygin");


    // Constructors

        //- Construct from mesh and thermo-mechanical properties
        MalyginMOXCreepModel
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );


    //- Destructor
    ~MalyginMOXCreepModel();


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