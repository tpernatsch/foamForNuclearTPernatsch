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
    Foam::LimbackCreepModel

Description
    Limback creep model

    \htmlinclude creepModel_LimbackCreepModel.html
    
SourceFiles
    LimbackCreepModel.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef LimbackCreepModel_H
#define LimbackCreepModel_H

#include "creepModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermoMechanicsSolver Declaration
\*---------------------------------------------------------------------------*/

class LimbackCreepModel
:
    public creepModel
{
    // Private Member Functions

        //- Name used for the fastFluence field
        const word fastFluenceName_;

        //- Name used for the fast flux field
        const word fastFluxName_;
        
        //- Pointer to fastFluence field
        const volScalarField* fastFluence_;
        
        //- Pointer to fast flux field
        const volScalarField* fastFlux_;

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

        //- Relaxation factor for creep strain
        scalar relax_;

        //- Type of cladding material for creep correlations. Default is SRA
        word cladType_;

        //- Gas Constant [ J/mol/K]
        const scalar R;
        
        //- Dimensionless constant for irradiation creep
        scalar C0; 
        
        //- Dimensionless constant for irradiation creep
        scalar C1;
        
        //- Dimensionless constant for irradiation creep
        scalar C2; 
        
        //- Material constant [K/Mpa/hr]
        scalar A;
        
        //- Activation energy per mol [J/mol]
        scalar Q;
        
        //- Primary creep exponent 
        scalar n;
        
        //- Dimensionless constants for thermal creep
        const scalar A1;
        
        //- Dimensionless constants for thermal creep
        const scalar a;
        
        //- Dimensionless constants for thermal creep
        const scalar A3;
        
        //- Material constant [n/cm2]^-A3
        const scalar A2;
        
        //- Dimensionless constants for primary creep
        const scalar C;
        
        //- Dimensionless constants for primary creep
        const scalar b;
        
        //- Dimensionless constants for primary creep
        const scalar d;
        
        //-Maerial constant [hr^b]
        const scalar B;
        
        //-Maerial constant [hr]
        const scalar D;
        
        //- Disallow default bitwise copy construct
        LimbackCreepModel(const LimbackCreepModel&);

        //- Disallow default bitwise assignment
        void operator=(const LimbackCreepModel&);

public:
 
    //- Runtime type information
    TypeName("LimbackCreepModel");

    // Constructors

        //- Construct from mesh and thermo-mechanical properties
        LimbackCreepModel
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );

    //- Destructor
    ~LimbackCreepModel();


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
