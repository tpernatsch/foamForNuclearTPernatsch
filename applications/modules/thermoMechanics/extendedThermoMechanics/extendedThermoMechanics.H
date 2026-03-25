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

#include "globalOptions.H"
#include "offbeatTime.H"
#include "thermalSubSolver.H"
#include "mechanicsSubSolver.H"
#include "gapGasModel.H"
#include "fissionGasRelease.H"
#include "burnup.H"
#include "neutronicsSubSolver.H"
#include "heatSource.H"
#include "fastFlux.H"
#include "sliceMapper.H"
#include "materials.H"
#include "rheology.H"
#include "corrosion.H"
#include "elementTransport.H"


#ifdef OPENFOAMFOUNDATION
#include "flowSubSolver.H"
#include "helperFuns.H"
#endif

#include "InterpolateTables.H"
#include "solver.H"

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
    public solver
{
    // Private data

        //- Vustom time
        const offbeatTime& runTime_;

        //- Main dictionary
        IOdictionary solverDict_;
        
        //- Global options (e.g. pin direction)
        globalOptions globalOpt_;

        //- Materials
        autoPtr<materials> mat_;

        //- Slice mapper
        autoPtr<sliceMapper> mapper_;

        //- Rheology
        autoPtr<rheology> rheo_;

        //- Heat source
        autoPtr<heatSource> heatSrc_;

        //- Burnup and depletion
        autoPtr<burnup> burnup_;

        //- Neutronics solver
        autoPtr<neutronicsSubSolver> neutronics_;

        //- Fast flux
        autoPtr<fastFlux> fastFlux_;

        //- Thermal solver
        autoPtr<thermalSubSolver> thermal_;

        //- Mechanics solver
        autoPtr<mechanicsSubSolver> mechanics_;

        //- Fission gas release
        autoPtr<fissionGasRelease> fgr_;

        //- Gap gas model
        autoPtr<gapGasModel> gapGas_;

        //- Corrosion model
        autoPtr<corrosion> corrosion_;

        //- Element transport solvers
        autoPtr<elementTransport> elementTransport_;

#ifdef OPENFOAMFOUNDATION
        //- Flow Solver
        autoPtr<flowSubSolver> flow_;
#endif

        //- Number of Outer iterations (reset at each time steo)
        scalar nOuterIter_;

        //- Total number of iterations from start to finish
        scalar totalIterations_;

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
            dynamicFvMesh& mesh
        );

    // Selectors

        //- Destructor
        virtual ~extendedThermoMechanics()
        {}

    // Member Functions

        //- Run extendedThermoMechanics
        #ifdef isCommDataLayerIncluded
        virtual void correctPhysics(fmuControl* fmu);
        #else
        virtual void correctPhysics();
        #endif

        //- Return next deltaT according to adjustable time class
        virtual scalar maxDeltaT();

        //- Return end time
        virtual scalar endTime();

    // Access Functions
        
        //- Return total number of iterations
        virtual scalar totalIterations()
        {
            return totalIterations_;
        }
        
        //- Return current number of outer iterations
        virtual scalar nOuterIter()
        {
            return nOuterIter_;
        }

        //- Give back current residual
        // TODO: change to obtain real residual
        virtual scalar getResidual()
        {
            return 0;
        };

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} //End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
