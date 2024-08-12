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
    Foam::offbeat

Description
    

SourceFiles
    offbeat.C

\*---------------------------------------------------------------------------*/

#ifndef offbeat_H
#define offbeat_H

// Include standard classes
#include "fvCFD.H"
#include "SortableList.H"

// Include main offbeat classes
#include "globalFieldLists.H"
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
#include "timeHandling.H"
#include "sliceMapper.H"
#include "materials.H"
#include "rheology.H"
#include "corrosion.H"
#include "userParameters.H"
#include "elementTransport.H"

#include "InterpolateTables.H"
#include "solver.H"

#include "runTimeSelectionTables.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{


/*---------------------------------------------------------------------------*\
                         Class offbeat Declaration
\*---------------------------------------------------------------------------*/

class offbeat
:
    public solver
{
    // Private data

        const offbeatTime& runTime_;
        fvMesh& mesh_;
        IOdictionary solverDict_;
        timeHandling adjustableTime_;
        globalOptions globalOpt_;
        autoPtr<materials> mat_;
        autoPtr<sliceMapper> mapper_;
        autoPtr<rheology> rheo_;
        autoPtr<heatSource> heatSrc_;
        autoPtr<burnup> burnup_;
        autoPtr<neutronicsSubSolver> neutronics_;
        autoPtr<fastFlux> fastFlux_;
        autoPtr<thermalSubSolver> thermal_;
        autoPtr<mechanicsSubSolver> mechanics_;
        autoPtr<fissionGasRelease> fgr_;
        autoPtr<gapGasModel> gapGas_;
        autoPtr<corrosion> corrosion_;
        autoPtr<elementTransport> elementTransport_;
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        offbeat(const offbeat&);

        //- Disallow default bitwise assignment
        void operator=(const offbeat&);

public:

    TypeName("offbeat");
    
    // Constructors

    

        //- Construct from mesh, materials and dictionary
        offbeat
        (
            fvMesh& mesh
        );

    // Selectors

        //- Destructor
        virtual ~offbeat()
        {}

    // Member Functions

        // Run offbeat
        virtual void correctPhysics();

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


        virtual autoPtr<fvMesh> createBaffleLessMesh();
        // virtual PtrList<volScalarField>  createScalarBaffleLessFields(){};
        // virtual PtrList<volVectorField>  createVectorBaffleLessFields(){};
        virtual void correctBaffleLessFields(){};

        virtual void deformMesh(){};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} //End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
