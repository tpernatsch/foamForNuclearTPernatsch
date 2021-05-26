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
    This class encapsulates a two-fluid phases + 1 stationary structure 
    Eulerian solver with a user-defineable regime-map based phyiscs modelling.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    twoPhase.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhase_H
#define twoPhase_H

#include "thermalHydraulicModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class regime;
class regimeMapModel;
class phaseChangeModel;

namespace thermalHydraulicModels
{

/*---------------------------------------------------------------------------*\
                      Class twoPhase Declaration
\*---------------------------------------------------------------------------*/

class twoPhase
:
    public thermalHydraulicModel
{
public :

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
        autoPtr<structureModel> structurePtr_;
        structureModel& structure_;

        //- Fluid1
        fluid fluid1_;

        //- Fluid2
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

        //- Reynolds of the two-phase mixture. Logically this belongs to none
        //  of the pairs above, as it depends on both fluids and the structure 
        //  (namely its hydraluic diameter) at the same time
        volScalarField ReTwoPhase_;

    //- Mixture quantities (while the solver is Euler-Euler, these can be of
    //  intereset and are occasionally used)
        
        //- Mixture velocity
        volVectorField U_;
        
        //- Mixture density
        volScalarField rho_;
    
    //- Phase interaction related  
        
        //- table of drag coefficients between:
        //      - fluid1 and fluid2;
        //      - fluid1 and the structure
        //      - fluid2 and the structure
        volTensorFieldPtrTable Kds_;

        //- table of heat transfer coefficients between:
        //      - fluid1 and interface between fluid1 and fluid2;
        //      - fluid2 and interface between fluid2 and fluid1;
        //      - fluid1 and the structure
        //      - fluid2 and the structure
        volScalarFieldPtrTable htcs_;

        //- Virtual mass coefficient, currently handled outside the
        //  regime framework
        scalar Vm_;

        //- interfacial area between fluid1 and fluid2
        volScalarField iA12_;

        //- Mass transfer term from fluid1 to fluid2 so that:
        /*
            ddt(rho1*alpha1) + div(alphaRhoPhi1) = -dmdt
            ddt(rho2*alpha2) + div(alphaRhoPhi2) = dmdt
        */
        volScalarField dmdt_;

        //- Interfacial temperature between phase1 and phase2
        volScalarField iT12_;

        //- Regime map model
        autoPtr<regimeMapModel> regimeMap_;

        //- Phase change model
        autoPtr<phaseChangeModel> phaseChange_;

    //- Extra field for under-relaxation of certain quantities
        fvVectorMatrix  VmVolForce1_;
        fvVectorMatrix  VmVolForce2_;

    //- Miscellanea

    //- Global flags for time/subCycle control

        //- Self explanatory, come on...
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
            myPimpleControl& pimple, 
            fv::options& fvOptions
        );

        //- Destructor
        virtual ~twoPhase(){}

    //- Protected member Functions

protected :

    //- Under-relax coupling fields (i.e. Kd, htc) in the provided
    //  table
    template<class fieldType>
    void relaxFieldPtrTable
    (
        const word& relaxFactorName,
        HashPtrTable
        <
            fieldType, 
            word, 
            word::hash
        >& fieldPtrTable
    );
    
    //- Public member Functions

public :

        //- Return the mixture density
        tmp<volScalarField> rho() const;

        //- Return the mixture velocity
        tmp<volVectorField> U() const;

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

        //- Correct alphas, pressure, velocities without updating 
        //  the flow regime map
        virtual void correctFluidMechanics(scalar& residual);

        //- Correct flow regime map and all the coupling fields
        virtual void correctRegimes
        (
            bool solveFluidDynamics, 
            bool solveEnergy
        );

        //- Update and print Courant number
        virtual void correctCourant();

        //- Custom adjustTimeStep
        virtual void adjustTimeStep();

        //- Update continuity errors of the fluids
        void correctContErrs();

        //- Print continuity error information
        void printContErrs();

        //- Computes the cumulated continuity errors
        void calcCumulContErrs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalHydraulicModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
