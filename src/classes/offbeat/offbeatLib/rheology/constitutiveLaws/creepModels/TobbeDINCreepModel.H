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
    Foam::TobbeDINCreepModel

Description
    Creep model from Tobbe (1975) correlation, for the specific steel DIN 1.4970m of the 15-15 Ti class.
    \htmlinclude creepModel_TobbeDINCreepModel.html

SourceFiles
    TobbeDINCreepModel.C

\todo
    check that coding for boundary conditions is correctA. Scolaro

\*---------------------------------------------------------------------------*/

#ifndef TobbeDINCreepModel_H
#define TobbeDINCreepModel_H

#include "creepModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermoMechanicsSolver Declaration
\*---------------------------------------------------------------------------*/

class TobbeDINCreepModel
:
    public creepModel
{
    // Private Member Functions

        //- Name used for the fast flux field
        const word fastFluxName_;
        
        //- Pointer to fast flux field
        const volScalarField* fastFlux_;

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
        
        scalar timeIndex_;

        scalar relax_;

        //- Constants
        const scalar par1;
        const scalar par2;
        const scalar par3;
        const scalar par4;
        const scalar par5;
        
        //- Disallow default bitwise copy construct
        TobbeDINCreepModel(const TobbeDINCreepModel&);

        //- Disallow default bitwise assignment
        void operator=(const TobbeDINCreepModel&);

public:
 
    //- Runtime type information
    TypeName("TobbeDINCreepModel");

    // Constructors

        //- Construct from mesh and thermo-mechanical properties
        TobbeDINCreepModel
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );

    //- Destructor
    ~TobbeDINCreepModel();


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
