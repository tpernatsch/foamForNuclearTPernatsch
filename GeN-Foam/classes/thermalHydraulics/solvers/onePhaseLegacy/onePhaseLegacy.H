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
    Foam::onePhaseLegacy

Description
    This class encapsulates a single fluid + 1 stationary structure 
    Eulerian solver with a user-defineable regime-map based phyiscs modelling.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    onePhaseLegacy.C

\*---------------------------------------------------------------------------*/

#ifndef onePhaseLegacy_H
#define onePhaseLegacy_H

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
                      Class onePhaseLegacy Declaration
\*---------------------------------------------------------------------------*/

class onePhaseLegacy
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

    //- Additional, temporary fields for the solver

        //- Darcy velocity which is the 
        autoPtr<volVectorField> UDarcy_;

        //- Continuous Darcy flux 
        autoPtr<surfaceScalarField> phiDarcy_;

        //- The temperature correction term for the Boussinesq approximation
        autoPtr<volScalarField> rhok_;

    //- Additional dimensioned constants for the Boussinesq approximation

        //- Reference temperature
        autoPtr<dimensionedScalar> Tref_;

        //- Thermal expansion coefficient
        autoPtr<dimensionedScalar> beta_;

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
        volTensorFieldPtrTable Kds_;

        //- table of heat transfer coefficients between the fluid and each
        //  structure region
        volScalarFieldPtrTable htcs_;

        //- Regime map
        autoPtr<regimeMapModel> regimeMap_;

    //- Container for the cumulative continuity error
        scalar cumulativeContErr_;

public:

    TypeName("onePhaseLegacy");

    //- Constructors

        onePhaseLegacy
        (
            Time& time,
            fvMesh& mesh, 
            myPimpleControl& pimple, 
            fv::options& fvOptions
        );

        //- Destructor
        virtual ~onePhaseLegacy(){}

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
