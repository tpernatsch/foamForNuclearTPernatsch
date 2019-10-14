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
    Class to model a generic fluid phase. It mainly acts as a variable
    placeholder, as little functionality is deferred to the class itself.

Author
    Stefan Radman (stefanradman92@gmail.com, stefan.radman@epfl.ch)

Class
    Foam::fluid

SourceFiles
    fluid.C

\*---------------------------------------------------------------------------*/

#ifndef fluid_H
#define fluid_H

#include "phaseBase.H"
#include "rhoThermo.H"
#include "fluidThermo.H"
#include "PhaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
template<class fluidThermo>
class PhaseCompressibleTurbulenceModel;

/*---------------------------------------------------------------------------*\
                           Class fluid Declaration
\*---------------------------------------------------------------------------*/

class fluid
:
    public phaseBase
{
private :

    //- Thermophysical properties
    autoPtr<rhoThermo> thermo_;

    //- Velocity
    volVectorField U_;

    //- Velocity magnitude;
    volScalarField magU_;

    //- Superficial flux of the phase
    autoPtr<surfaceScalarField> phiPtr_;

    //- Volumetric flux of the phase
    surfaceScalarField alphaPhi_;

    //- Mass flux of the phase
    surfaceScalarField alphaRhoPhi_;

    //- Dilation term of the phase.
    //  It encompasses the "compressible" part of the continuity
    //  equation, including mass sources/sinks, if any.
    //  Example given in meta-code:
    /*
        ddt(rho*alpha) + div(alphaRhoPhi) = dmdt

        -> (expand the derivatives)

        rho1*ddt(alpha) + rho1*div(alphaPhi) +
        alpha1*ddt(rho) + alphaPhi1&grad(rho) = dmdt

        -> (re-arrange)

        ddt(alpha) + div(alphaPhi) 
        =
        (1/rho)*
        (
            alpha*ddt(rho) + alphaPhil&grad(rho) + dmdt
        ) 
        =
        dgdt
    */
    volScalarField dgdt_;

    //- Continuity error of the phase 
    //  (== fvc::ddt(alpha, rho) + fvc::div(alphaRhoPhi)) [kg/m3/s]
    volScalarField contErr_;

    //- Energy density error of the phase [J/m3]
    volScalarField heErr_;

    //- Characteristic dimension of the fluid. This is equal to the structure
    //  hydraulic diameter if this fluid is continuous, else it is specified
    //  by the fluidDiameterModel within fluidGeometry. Information on whether
    //  a fluid is continuous or dispersed in a certain regime is specified in
    //  the fluidGeometry dict as well. Since this is regime-dependent, this 
    //  field is manipulated by reference by each regimeModel and it is never
    //  directly modified by the fluid class itself
    volScalarField Dh_;

    //- Marker field that is 1 where the phase is dispersed, 0 where it is 
    //  continuous, and between 0 and 1 in interpolation regimes
    scalarField dispersion_;

    //- Turbulence model
    typedef PhaseCompressibleTurbulenceModel<fluidThermo> FFSTurbulenceModel;
    autoPtr<FFSTurbulenceModel> turbulence_;

    //- Phase fraction below which the fluid energy equation is not solved
    //  and the phase temperature is set to the interfacial temperature
    dimensionedScalar thermoResidualAlpha_;

public:

    // Constructors

        fluid
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& phaseName
        );

        //- Destructor
        virtual ~fluid();

    // Member Functions

        //- To remove abstractness from phaseBase
        virtual void abstract() {}

        //- Ensure that the flux at inflow/outflow BCs is preserved
        void correctInflowOutflow(surfaceScalarField& alphaPhi) const;

        // Access

            const FFSTurbulenceModel& turbulence() const
            {
                return turbulence_();
            }

            const rhoThermo& thermo() const
            {
                return thermo_();
            }

            rhoThermo& thermo()
            {
                return thermo_();
            }

            const volScalarField& T() const
            {
                return thermo_->T();
            }

            volScalarField& T()
            {
                return thermo_->T();
            }

            const volScalarField& rho() const
            {
                return thermo_->rho();
            }

            const volScalarField& Cp() const
            {
                return thermo_->Cp();
            }

            const volScalarField& Cv() const
            {
                return thermo_->Cv();
            }

            const volVectorField& U() const
            {
                return U_;
            }

            volVectorField& U()
            {
                return U_;
            }

            const volScalarField& magU() const
            {
                return magU_;
            }

            volScalarField& magU()
            {
                return magU_;
            }

            const surfaceScalarField& phi() const
            {
                return phiPtr_();
            }

            surfaceScalarField& phi()
            {
                return phiPtr_();
            }

            const surfaceScalarField& alphaPhi() const
            {
                return alphaPhi_;
            }

            surfaceScalarField& alphaPhi()
            {
                return alphaPhi_;
            }

            const surfaceScalarField& alphaRhoPhi() const
            {
                return alphaRhoPhi_;
            }

            surfaceScalarField& alphaRhoPhi()
            {
                return alphaRhoPhi_;
            }

            const volScalarField& dgdt() const
            {
                return dgdt_;
            }

            volScalarField& dgdt()
            {
                return dgdt_;
            }

            const volScalarField& contErr() const
            {
                return contErr_;
            }

            volScalarField& contErr()
            {
                return contErr_;
            }

            const volScalarField& heErr() const
            {
                return heErr_;
            }

            volScalarField& heErr()
            {
                return heErr_;
            }

            const volScalarField& Dh() const
            {
                return Dh_;
            }

            volScalarField& Dh()
            {
                return Dh_;
            }

            const scalarField& dispersion() const
            {
                return dispersion_;
            }

            scalarField& dispersion()
            {
                return dispersion_;
            }

            const dimensionedScalar thermoResidualAlpha() const
            {
                return thermoResidualAlpha_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
