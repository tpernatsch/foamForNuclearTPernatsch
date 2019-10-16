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
class saturationModel;

namespace thermalHydraulicModels
{

/*---------------------------------------------------------------------------*\
                      Class twoPhase Declaration
\*---------------------------------------------------------------------------*/

class twoPhase
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

        //- Fluid1
        fluid fluid1_;

        //- Fluid2
        fluid fluid2_;

        //- Structure
        autoPtr<structureModel> structurePtr_;
        structureModel& structure_;

        //- Total fluid phase fraction
        volScalarField movingAlpha_;

    //- Phase pairs
        
        //- Fluid1 - fluid2
        FFPair FFPair_;
        
        //- Fluid1 - structure
        FSPair F1SPair_;
        
        //- Fluid2 - structure
        FSPair F2SPair_;
    
    //- Phase interaction related  
        
        //- table of drag coefficients between:
        //      - fluid1 and fluid2;
        //      - fluid1 and the structure
        //      - fluid2 and the structure
        volTensorFieldTable Kds_;

        //- Avoid momentum decoupling between the two fluids by imposing a 
        //  minimum value for the drag coefficient
        dimensionedScalar residualKd_;

        //- table of heat transfer coefficients between:
        //      - fluid1 and interface between fluid1 and fluid2;
        //      - fluid2 and interface between fluid2 and fluid1;
        //      - fluid1 and the structure
        //      - fluid2 and the structure
        volScalarFieldTable htcs_;

        //- interfacial area between fluid1 and fluid2
        volScalarField iA12_;

        //- fraction of interfacialArea of the corresponding structure present
        //  in a given cell that is in contact with fluid1
        volScalarField frac1_;

        //- fraction of interfacialArea of the corresponding structure present
        //  in a given cell that is in contact with fluid2
        volScalarField frac2_;

        //- Mass transfer term from fluid1 to fluid2 so that:
        /*
            ddt(rho1*alpha1) + div(alphaRhoPhi1) = -dmdt
            ddt(rho2*alpha2) + div(alphaRhoPhi2) = dmdt
        */
        volScalarField dmdt_;

        //- Interfacial temperature between phase1 and phase2
        volScalarField iT12_;

        //- Saturation model
        autoPtr<saturationModel> saturation_;

        //- Regime map
        autoPtr<regimeMapModel> regimeMap_;

    //- Solution control

        //- Enum for selected phaseChange model
        const phaseChangeModel phaseChangeModel_;
        const bool phaseChange_;
        const word volatilePhaseName_;

        //- Partial elimination flag
        const Switch partialElimination_;

        //- Face momentum flag
        const Switch faceMomentum_;

public:

    TypeName("twoPhase");

    //- Constructors

        twoPhase
        (
            const fvMesh& mesh, 
            pimpleControl& pimple, 
            fv::options& fvOptions
        );

        //- Destructor
        virtual ~twoPhase(){}

    //- Member Functions

        //- Return the mixture density
        tmp<volScalarField> rho() const;

        //- Return the mixture velocity
        tmp<volVectorField> U() const;

        //- Solve all physics altogether
        virtual void correct(scalar& residual);

        //- Correct fluid and structure temperatures without
        //  updating the flow regime map
        virtual void correctEnergy(scalar& residual);

        //- Correct alphas, pressure, velocities without updating 
        //  the flow regime map
        virtual void correctFluidMechanics(scalar& residual);

        //- Correct flow regime map and all the coupling fields
        virtual void correctRegimes();

        //- Update and print Courant number
        virtual void correctCourant();

        //- Update continuity errors of the fluids
        void correctContErrs();

        //- GeN-Foam compatibility

            virtual const volScalarField& alphaLiquid() const
            {
                return fluid1_;
            }

            virtual const volScalarField& alphaVapour() const
            {
                return fluid2_;
            }

            virtual const volScalarField& alphaStructure() const
            {
                return structure_;
            }

            virtual const volScalarField& TLiquid() const
            {
                return fluid1_.thermo().T();
            }

            virtual const volScalarField& TVapour() const
            {
                return fluid2_.thermo().T();
            }

            virtual const volScalarField& TStructure() const
            {
                return structure_.Tpas();
            }

            virtual const volVectorField& ULiquid() const
            {
                return fluid1_.U();
            }

            virtual const volVectorField& UVapour() const
            {
                return fluid2_.U();
            }

            virtual tmp<volScalarField> rhoLiquid() const
            {
                return fluid1_.thermo().rho();
            }

            virtual tmp<volScalarField> rhoVapour() const
            {
                return fluid2_.thermo().rho();
            }

            virtual const volScalarField& TavFuel() const
            {
                return structure_.TavFuel();
            }

            virtual const volScalarField& TavClad() const
            {
                return structure_.TavClad();
            }

            virtual volScalarField& nuclearFuelPowerDensity()
            {
                return structure_.nuclearFuelPowerDensity();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalHydraulicModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
