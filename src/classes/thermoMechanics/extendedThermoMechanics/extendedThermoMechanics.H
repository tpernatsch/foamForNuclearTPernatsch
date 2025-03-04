/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::extendedThermoMechanics

Description


SourceFiles
    extendedThermoMechanics.C

\*---------------------------------------------------------------------------*/

#ifndef extendedThermoMechanics_H
#define extendedThermoMechanics_H

// Include standard classes
#include "fvCFD.H"
#include "SortableList.H"

#include "globalFieldLists.H"
#include "globalOptions.H"
#include "offbeatTime.H"
#include "thermalSubSolver.H"
#include "mechanicsSubSolver.H"
#include "heatSource.H"
#include "timeHandling.H"
#include "sliceMapper.H"
#include "materials.H"
#include "rheology.H"
#include "userParameters.H"

#include "InterpolateTables.H"
#include "thermoMechanics.H"

#include "runTimeSelectionTables.H"
#include "dynamicFvMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{


/*---------------------------------------------------------------------------*\
                         Class extendedThermoMechanics Declaration
\*---------------------------------------------------------------------------*/

class extendedThermoMechanics
:
    public thermoMechanics
{
    // Private data

        const offbeatTime& runTime_;
        dynamicFvMesh& mesh_;
        IOdictionary solverDict_;
        timeHandling adjustableTime_;
        globalOptions globalOpt_;
        autoPtr<materials> mat_;
        autoPtr<sliceMapper> mapper_;
        autoPtr<rheology> rheo_;
        autoPtr<heatSource> heatSrc_;
        autoPtr<thermalSubSolver> thermal_;
        autoPtr<mechanicsSubSolver> mechanics_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        extendedThermoMechanics(const extendedThermoMechanics&);

        //- Disallow default bitwise assignment
        void operator=(const extendedThermoMechanics&);

public:

    TypeName("extendedThermoMechanics");

    // Constructors



        //- Construct from mesh, materials and dictionary
        extendedThermoMechanics
        (
            fvMesh& mesh
        );

    // Selectors

        //- Destructor
        virtual ~extendedThermoMechanics()
        {}

    // Member Functions

        // Run extendedThermoMechanics
        virtual void correctPhysics();
        virtual void correctTightlyCoupledPhysics();

        // Set next deltaT according to adjustable time class
        virtual scalar maxDeltaT();


        // Returns ref to runTime
        inline virtual const Time& getRunTime()
        {
            return runTime_;
        }

        // Returns max value of T at the outer cladding patch
        virtual scalar getMaxToutClad
        (
            const word& cladOuterPatchName
        );

        inline virtual const fvMesh& mesh()
        {
            return mesh_;
        }

        virtual scalar getResidual(){return 0;};


        virtual void correctBaffleLessFields(){};

        virtual void deformMesh(){};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} //End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
