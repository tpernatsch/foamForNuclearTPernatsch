/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::thermalHydraulicModel

Description
    Run-time selectable base abstract class to represent a thermal-hydraulic 
    solver for an indefinite number of phases.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    thermalHydraulicModel.C

\*---------------------------------------------------------------------------*/

#ifndef thermalHydraulicModel_H
#define thermalHydraulicModel_H

#include "fvCFD.H"
#include "fvOptions.H"
#include "IOdictionary.H"
#include "pimpleControl.H"
#include "runTimeSelectionTables.H"

//- It is pointless to include these in the .C files as any derived solver 
//  class will need all of these anyway, so I include them in the header too
#include "volFields.H"
#include "UniformField.H"
#include "surfaceFields.H"
#include "surfaceInterpolate.H"
#include "fixedValueFvsPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvmDdt.H"
#include "fvcCurl.H"
#include "fvcFlux.H"
#include "fvMatrix.H"
#include "fvcSnGrad.H"
#include "fvmLaplacian.H"
#include "HashPtrTable.H"

#include "fluid.H"
#include "structureModel.H"
#include "FFPair.H"
#include "FSPair.H"
#include "typedefs.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class thermalHydraulicModel Declaration
\*---------------------------------------------------------------------------*/

class thermalHydraulicModel
:
    public IOdictionary
{
protected:

    Time& runTime_;

    fvMesh& mesh_;

    pimpleControl& pimple_;

    fv::options& fvOptions_;

    //- Multiple Rotating reference Frames
    IOMRFZoneList MRF_;

    //- Total volumetric flux
    surfaceScalarField phi_;

    //- Gravity related fields
    uniformDimensionedVectorField g_;
    uniformDimensionedScalarField hRef_;
    volScalarField gh_;
    surfaceScalarField ghf_;

    //- Pressure
    volScalarField p_rgh_;
    volScalarField p_;

    //- Minimum allowable pressure
    dimensionedScalar pMin_;

    //- Reference cell and pressure if no fixedValue BCs provided
    label pRefCell_;
    scalar pRefValue_;

    //- Pressure power term
    volScalarField dpdt_;

    //- Initial fluid(s) mass
    dimensionedScalar initialFluidMass_;

    //- Maximum Courant number
    scalar CoNum_;

    //- Mean Courant number
    scalar meanCoNum_;

public:

    //- Runtime type information
    TypeName("thermalHydraulicModel");

    //- Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            thermalHydraulicModel,
            thermalHydraulicModels,
            (
                Time& time,
                fvMesh& mesh,
                pimpleControl& pimple,
                fv::options& fvOpt
            ),
            (time, mesh, pimple, fvOpt)
        );

    //- Constructors

        thermalHydraulicModel
        (
            Time& time,
            fvMesh& mesh,
            pimpleControl& pimple,
            fv::options& fvOpt
        );

        //- Destructor
        virtual ~thermalHydraulicModel();

    //- Selectors

        static autoPtr<thermalHydraulicModel> New
        (
            Time& time,
            fvMesh& mesh,
            pimpleControl& pimple,
            fv::options& fvOpt
        );

    //- Member Functions

        //- Solve thermal-hydraulics, combination of the other corrects
        virtual void correct
        (
            scalar& residual,
            bool solveFluidDynamics, 
            bool solveEnergy
        ) = 0;

        //- Correct fluid and structure temperatures without
        //  updating the flow regime map
        virtual void correctEnergy(scalar& residual) = 0;

        //- Correct alphas, pressure, velocities without updating 
        //  the flow regime map
        virtual void correctFluidMechanics(scalar& residual) = 0;

        //- Correct flow regime map and all the coupling fields
        virtual void correctRegimes
        (
            bool solveFluidDynamics, 
            bool solveEnergy
        ) = 0;

        //- Update and print Courant number
        virtual void correctCourant() = 0;

        //- Stop execution (for development purposes)
        void stop();

        //- A function to allow for the implementation of custom adjustTimeStep
        //  criteria
        virtual void adjustTimeStep();

        // Access

            //- Return maximum Courant number
            scalar CoNum() const
            {
                return CoNum_;
            }

            //- Return mean Courant number
            scalar meanCoNum() const
            {
                return meanCoNum_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
