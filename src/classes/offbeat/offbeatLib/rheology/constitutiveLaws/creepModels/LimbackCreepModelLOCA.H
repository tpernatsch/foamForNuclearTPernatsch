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
    Foam::LimbackCreepModelLOCA 

Description
    Limback creep model extended to LOCA.

    \htmlinclude creepModel_LimbackCreepModelLOCA.html

SourceFiles 
    LimbackCreepModelLOCA.C

\mainauthor
    E. Brunetto, A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef LimbackCreepModelLOCA_H
#define LimbackCreepModelLOCA_H

#include "creepModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam 
{

/*---------------------------------------------------------------------------*\
                    Class "LimbackCreepModelLOCA" Declaration
\*---------------------------------------------------------------------------*/

class LimbackCreepModelLOCA
:
    public creepModel
{
    // Private Data

        //- Name used for the fastFluence field
        const word fastFluenceName_;
 
        //- Name used for the fast flux field
        const word fastFluxName_;
        
        //- Pointer to fastFluence field
        const volScalarField* fastFluence_;
        
        //- Pointer to fast flux field
        const volScalarField* fastFlux_;

        //- Pointer to beta fraction field
        const volScalarField* betaFraction_;

        //- Irradiation creep Strain 
        volSymmTensorField& epsilonCreepIrr_;

        //- Thermal creep Strain 
        volSymmTensorField& epsilonCreepTh_;

        //- Primary Creep Strain 
        volSymmTensorField& epsilonCreepPrim_;
        
        //- equivalent Irradiation Creep Strain 
        volScalarField& epsilonCreepIrrEq_;
        
        //- Increment of equivalent Irradiation Creep Strain 
        volScalarField& DepsilonCreepIrrEq_;
        
        //- Equivalent Secondary Thermal Creep Strain 
        volScalarField& epsilonCreepThEq_;
        
        //- Increment of equivalent Secondary Thermal Creep Strain 
        volScalarField& DepsilonCreepThEq_;
        
        //- Equivalent Primary Creep Strain 
        volScalarField& epsilonCreepPrimEq_;
        
        //- Increment of equivalent Primary Creep Strain 
        volScalarField& DepsilonCreepPrimEq_;

        scalar timeIndex_;

        scalar relax_;

        // Activate irradiation creep
        Switch irradiationCreep_;

        // Activate primary creep
        Switch primaryCreep_;

        // Activate Newton-Raphson method for convergence
        Switch NewtonRaphsonMethod_;

        // Tolerance for inner convergence of the Newton-Raphson method
        scalar NewtonRaphsonTolerance_;

        //- Dimensionless constant
        const scalar C0; 
        const scalar C1;
        const scalar C2; 
        
    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        LimbackCreepModelLOCA(const LimbackCreepModelLOCA&);

        //- Disallow default bitwise assignment
        void operator=(const LimbackCreepModelLOCA&);

public:
 
    //- Runtime type information
    TypeName("LimbackCreepModelLOCA");


    // Constructors

        //- Construct from mesh and thermo-mechanical properties
        LimbackCreepModelLOCA
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );


    //- Destructor
    ~LimbackCreepModelLOCA();


    // Member Functions

        //- Update mechanical properties according to rheology model        
        virtual void correctCreep
        (
            volSymmTensorField& epsilonEl,
            const labelList& addr
        );

        //- Compute irradiation creep rate
        virtual scalar calcIrradiationCreepRate
        (
            const scalar& flux,
            const scalar& sigmaEff
        );

        //- Compute derivative of irradiation creep rate
        virtual scalar calcIrradiationCreepRatePrime
        (
            const scalar& flux,
            const scalar& sigmaEff,
            const scalar& creepIrrRate
        );

        //- Compute secondary creep rate in normal operating conditions
        virtual Foam::scalar calcNormalThermalCreepRate
        (
            const scalar& fastFluence, 
            const scalar& T, 
            const scalar& sigmaEff
        );

        //- Compute derivative of secondary creep rate in normal operating conditions
        virtual Foam::scalar calcNormalThermalCreepRatePrime
        (
            const scalar& fastFluence, 
            const scalar& T,
            const scalar& sigmaEff, 
            const scalar& creepThRate
        );

        //- Compute primary creep rate
        virtual Foam::scalar calcPrimaryCreepRate
        (
            const scalar& creepIrrRate,
            const scalar& creepThRate,
            const scalar& epsilonCreepPrimEqOld
        );

        //- Compute derivative of primary creep rate
        virtual Foam::scalar calcPrimaryCreepRatePrime
        (
            const scalar& creepIrrRate,
            const scalar& creepThRate,
            const scalar& creepIrrRatePrime,
            const scalar& creepThRatePrime,
            const scalar& epsilonCreepPrimEq,
            const scalar& epsilonCreepPrimEqOld
        );
        
        //- Compute secondary creep rate in LOCA conditions
        virtual Foam::scalar calcLOCAThermalCreepRate
        (
            const scalar& T, 
            const scalar& betaFrac, 
            const scalar& sigmaEff, 
            scalar creepThRateOld
        );
        
        //- Compute derivative of secondary creep rate in LOCA conditions
        virtual Foam::scalar calcLOCAThermalCreepRatePrime
        (
            const scalar& T, 
            const scalar& betaFrac,
            const scalar& sigmaEff,
            const scalar& creepRate
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //