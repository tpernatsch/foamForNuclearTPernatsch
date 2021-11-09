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
    Foam::twoPhase

Description
    This class encapsulates a two fluid + 1 stationary structure 
    Eulerian solver

Author
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch / )

SourceFiles
    twoPhase.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhase_H
#define twoPhase_H

#include "thermalHydraulicsModel.H"
#include "twoPhaseDragFactor.H"
#include "myOps.H"
#include "upwind.H"
#include "gaussConvectionScheme.H"
#include "myMULES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace thermalHydraulicsModels
{

/*---------------------------------------------------------------------------*\
                      Class twoPhase Declaration
\*---------------------------------------------------------------------------*/

class twoPhase
:
    public thermalHydraulicsModel
{
public:

    enum class alphaEqnsSolver
    {
        MULES,
        implicitUpwind
    };

    static const Enum<alphaEqnsSolver> alphaEqnsSolverNames_;

    enum class partialEliminationMode
    {
        none, 
        legacy,
        implicit,
        implicitWithDmdt
    };

    static const Enum<partialEliminationMode> partialEliminationModeNames_;

    enum class contErrCompensationMode
    {
        none,
        Su,
        Sp, 
        SuSp
    };

    static const Enum<contErrCompensationMode> contErrCompensationModeNames_;

    enum class heStabilizationMode
    {
        cutoff,
        source
    };

    static const Enum<heStabilizationMode> heStabilizationModeNames_;

private:

    //- Phases

        //- Structure
        structure structure_;

        //- Fluids
        fluid fluid1_;
        fluid fluid2_;

        //- Total fluid phase fraction
        volScalarField movingAlpha_;

    //- Phase pairs
        
        //- Fluid1 - fluid2
        FFPair FFPair_;
        
        //- Fluid1 - structure
        FSPair F1SPair_;
        
        //- Fluid2 - structure
        FSPair F2SPair_;

    //- 
        autoPtr<twoPhaseDragFactor> twoPhaseDragFactorPtr_;

    //- Mixture quantities

        //- Mixture velocity
        volVectorField U_;
        
        //- Mixture density
        volScalarField rho_;

    //- Miscellanea

    //- Global flags for time/subCycle control

        //- Self explanatory
        bool bothPhasesArePresent_;

        //- True if the temperature of of one of the two fluids is within a 
        //  certain margin from saturation temperature
        bool withinMarginToPhaseChange_;

    //- Solution control

        //-
        const alphaEqnsSolver alphaEqnsSolver_;

        //- partialElimination mode
        const partialEliminationMode partialEliminationMode_;

        //- contErrCompensation mode
        const contErrCompensationMode UContErrCompensationMode_;
        const contErrCompensationMode heContErrCompensationMode_;

        //-
        const heStabilizationMode heStabilizationMode_;

        //- Feature that removes velocity oscillations and/or smearing in the 
        //  main flow direction. However, if dealing with structures that have
        //  large porosity gradiens (i.e. > 0.5 between consecutive cells) or
        //  large gradients in the structure drag tensor, this feature may
        //  actually induce velocity oscillations in the directions 
        //  perpendicular to the main flow direction (nedless to say, not in 
        //  1D). Only applicable to the face-momentum based algorithm.
        //  This factor weights the "traditional" faceMomentum velocity 
        //  reconstruction (the one which always causes smearing/oscillations)
        //  in the main flow direction) with a moving phase fraction-adjusted 
        //  velocity reconstruction (which may induced the aforementioned
        //  transversal oscillations). If set to 0, it is a pure "traditional"
        //  faceMomentum.
        const scalar oscillationLimiterFraction_;

public:

    TypeName("twoPhase");

    //- Constructors

        twoPhase
        (
            Time& time,
            fvMesh& mesh, 
            customPimpleControl& pimple, 
            fv::options& fvOptions
        );

        //- Destructor
        virtual ~twoPhase(){}

    //- Member Functions

        //- Solve thermal-hydraulics
        virtual void correct
        (
            scalar& residual,
            bool solveFluidDynamics, 
            bool solveEnergy
        );

        //- Correct fluid and structure temperatures without
        //  updating the flow regime map
        virtual void correctEnergy(scalar& residual);

        //- Correct pressure, velocities without updating 
        //  the flow regime map
        virtual void correctFluidMechanics(scalar& residual);

        //- Correct models and coupling fields
        virtual void correctModels
        (
            bool solveFluidDynamics, 
            bool solveEnergy
        );

        //- Update and print Courant number
        virtual void correctCourant();

        //- Custom adjustTimeStep
        virtual void adjustTimeStep();

        //- Update continuity errors of the fluid
        void correctContErrs();

        //- Print continuity error information
        void printContErrs();

        //- Computes the cumulated continuity error
        void calcCumulContErrs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalHydraulicsModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
