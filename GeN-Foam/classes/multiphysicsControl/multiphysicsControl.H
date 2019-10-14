/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Description
    Class to handle fvSolution dictionaries for the three regions.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

Class
    Foam::multiphysicsControl

SourceFiles
    multiphysicsControl.C

\*---------------------------------------------------------------------------*/

#ifndef multiphysicsControl_H
#define multiphysicsControl_H

#include "fvMesh.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class multiphysicsControl Declaration
\*---------------------------------------------------------------------------*/

class multiphysicsControl
: 
    public pimpleControl
{
protected :

    //- Ref to runTime
    const Time& runTime_;

    //- Ref to thermalHydraulics mesh
    fvMesh& thermalHydraulicMesh_;

    //- Ref to neutronics mesh
    fvMesh& neutronicMesh_;

    //- Ref to thermoMechanics mesh
    fvMesh& thermoMechanicMesh_;

    IOdictionary        topLevelDict_;
    const dictionary&   thermalHydraulicDict_;
    const dictionary&   neutronicDict_;
    const dictionary&   thermoMechanicDict_;

    //- Solution flags

    //- solveFluidDynamics exists in pimpleControl as solveFlow_;

    bool solveEnergy_;

    bool solveNeutronics_;

    bool solveThermoMechanics_;

    //- Global loop parameter (i.e. outer iterations)

        bool    tightlyCoupled_;
        scalar  timeStepResidual_;
        label   maxTimeStepIterations_;

    //- Variables related to thermalHydraulics are in the PIMPLE base

    //- Related to neutronics

        bool    integralPredictor_;
        bool    implicitPredictor_;
        bool    aitkenAcceleration_;
        scalar  neutronIterationResidual_;
        label   maxNeutronIterations_;

    //- Related to thermoMechanics

        bool    compactNormalStress_;
        label   nThermoMechanicCorrs_;
        label   thermoMechanicCorr_;
        scalar  D_;

public:

    TypeName("multiphysicsControl");

    //- Constructors

        multiphysicsControl
        (
            const Time& runTime,
            fvMesh& THMesh,
            fvMesh& NMesh,
            fvMesh& TMMesh
        );

        //- Destructor
        virtual ~multiphysicsControl() {}

    //- Member Functions

        //- Override PIMPLE read to cope with different file structure
        virtual void read();

        //- Access

            //- Override main dict path
            virtual const dictionary& dict() const
            {
                return topLevelDict_;
            }

            bool solveFluidDynamics() const
            {
                return solveFlow_;
            }

            bool solveEnergy() const
            {
                return solveEnergy_;
            }

            bool solveNeutronics() const
            {
                return solveNeutronics_;
            }

            bool solveThermoMechanics() const
            {
                return solveThermoMechanics_;
            }

            label nOuterCorrectors() const
            {
                return nCorrPIMPLE_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
