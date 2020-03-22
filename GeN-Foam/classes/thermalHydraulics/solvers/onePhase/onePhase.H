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
    Foam::onePhase

Description
    This class encapsulates a single fluid + 1 stationary structure 
    Eulerian solver with a user-defineable regime-map based phyiscs modelling.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    onePhase.C

\*---------------------------------------------------------------------------*/

#ifndef onePhase_H
#define onePhase_H

#include "thermalHydraulicModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class regime;
class regimeMapModel;
class saturationModel;

namespace thermalHydraulicModels
{

/*---------------------------------------------------------------------------*\
                      Class onePhase Declaration
\*---------------------------------------------------------------------------*/

class onePhase
:
    public thermalHydraulicModel
{
public:

    enum class phaseChangeModel
    {
        none,
        heatConductionLimited, 
        individualSideDriven
    };

    static const Enum<phaseChangeModel> phaseChangeModelNames;

private:

    //- Phases

        //- Structure
        autoPtr<structureModel> structurePtr_;
        structureModel& structure_;

        //- Fluid
        fluid fluid_;

    //- Phase pairs
        
        //- Fluid2 - structure
        FSPair FSPair_;
    
    //- Phase interaction related  
        
        //- table of drag coefficients between the fluid and each
        //  structure region (currently only defined in terms of cellZone)
        volTensorFieldTable Kds_;

        //- table of heat transfer coefficients between the fluid and each
        //  structure region
        volScalarFieldTable htcs_;

        //- Regime map
        autoPtr<regimeMapModel> regimeMap_;

    //- Solution control

        //- Face momentum flag
        const Switch faceMomentum_;

        //- Feature that removes velocity oscillations and/or smearing in the 
        //  main flow direction. However, if dealing with structures that have
        //  large porosity gradiens (i.e. > 0.5 between consecutive cells) or
        //  large gradients in the structure drag tensor, this feature may
        //  actually induce velocity oscillations in the directions 
        //  perpendicular to the main flow directions (nedless to say, not in 
        //  1D). Only applicable to the face-momentum based algorithm.
        //  This factor weights the "traditional" faceMomentum velocity 
        //  reconstruction (the one which always causes smearing/oscillations)
        //  in the main flow direction) with a moving phase fraction-adjusted 
        //  velocity reconstruction (which may induced the aforementioned
        //  transversal oscillations). If set to 0, it is a pure "traditional"
        //  faceMomentum.
        const scalar oscillationLimiterFraction_;

public:

    TypeName("onePhase");

    //- Constructors

        onePhase
        (
            Time& time,
            fvMesh& mesh, 
            pimpleControl& pimple, 
            fv::options& fvOptions
        );

        //- Destructor
        virtual ~onePhase(){}

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

        //- Correct flow regime map and all the coupling fields
        virtual void correctRegimes
        (
            bool solveFluidDynamics, 
            bool solveEnergy
        );

        //- Update and print Courant number
        virtual void correctCourant();

        //- Update continuity error of the fluid
        void correctContErr();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalHydraulicModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
