/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::solvers::msfrParcelFoam

Description
    Transient solver for buoyant, turbulent, incompressible flow (Boussinesq
    approximation) for MSFR applications, with Lagrangian tracking of
    nuclear (fission product) particles undergoing decay, decay heat
    release and thermophoresis, and a surface film model for particle
    deposition -- the complete, actively-developed MSFR solver.

    Unlike OldParcelFoam (the original conversion, kept as a validated
    reference/comparison implementation), neutronics is factored out
    entirely into a separate neutroRegion running foamForNuclear's own
    native neutronics module (applications/modules/neutronics, e.g.
    diffusionNeutronics) rather than being solved inline on this region's
    own mesh. Coupling with the neutronics region is entirely through
    system/regionsDict's mappings: this region receives the volumetric
    fission power density as the mapped fields
    powerDensityLiquid_/powerDensityStructure_ (matching the exact field
    names applications/modules/thermalHydraulics's own onePhase solver
    reads), and sends T back out for the neutronics region's cross-section
    feedback.

    Verified against tutorials/reactorCases/2D_cavityMSFR_Reactor{,_split}
    and the CNRS multi-physics MSFR benchmark (Tiberga et al. 2020).

SourceFiles
    msfrParcelFoam.C

\*---------------------------------------------------------------------------*/

#ifndef msfrParcelFoam_H
#define msfrParcelFoam_H

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "turbulentTransportModel.H"
#include "singlePhaseTransportModel.H"
#include "surfaceFilmModel.H"
#include "basicNuclearCloud.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "radiationModel.H"
#include "IOMRFZoneList.H"
#include "solver.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace solvers
{

/*---------------------------------------------------------------------------*\
                      Class msfrParcelFoam Declaration
\*---------------------------------------------------------------------------*/

class msfrParcelFoam
:
    public solver
{
private:

    pimpleControl pimple_;

    //- Flow fields

    volVectorField U_;
    surfaceScalarField phi_;
    autoPtr<surfaceVectorField> Uf_;

    singlePhaseTransportModel laminarTransport_;

    //- Boussinesq / thermal properties

    dimensionedScalar beta_;
    dimensionedScalar CpRef_;
    dimensionedScalar TRef_;
    dimensionedScalar Pr_;
    dimensionedScalar Prt_;

    volScalarField T_;
    volScalarField p_rgh_;
    volVectorField vort_;

    autoPtr<incompressible::turbulenceModel> turbulence_;

    dimensionedScalar rhoInfValue_;
    volScalarField rhoInf_;
    volScalarField muc_;
    volScalarField rhok_;
    volScalarField alphat_;

    //- Gravity related fields

    uniformDimensionedVectorField g_;
    uniformDimensionedScalarField hRef_;
    volScalarField gh_;
    surfaceScalarField ghf_;

    volScalarField p_;
    label pRefCell_;
    scalar pRefValue_;

    volScalarField Qdot_;

    IOMRFZoneList MRF_;

    autoPtr<radiation::radiationModel> radiation_;
    dimensionedScalar rhoCpRef_;

    fv::options& fvOptions_;

    //- Fission product cloud
    basicNuclearCloud FP_;

    //- Surface film
    autoPtr<regionModels::surfaceFilmModel> tsurfaceFilm_;
    regionModels::surfaceFilmModel& surfaceFilm_;

    bool solvePrimaryRegion_;

    //- Volumetric fission power density mapped in from the neutronics
    //  region every timestep (system/regionsDict's mappings block); zero
    //  until the first mapping call. Same field names as
    //  applications/modules/thermalHydraulics's onePhase solver reads, so
    //  the same neutroRegion mapping block works unmodified.
    volScalarField powerDensityLiquid_;
    volScalarField powerDensityStructure_;

    scalar cumulativeContErr_;


public:

    TypeName("msfrParcelFoam");

    //- Constructors

        msfrParcelFoam
        (
            dynamicFvMesh& mesh
        );

        //- Destructor
        virtual ~msfrParcelFoam(){}

    //- Member Functions

        //- Solve thermal-hydraulics
        virtual void correctPhysics();

        virtual scalar maxDeltaT();

        virtual scalar getResidual()
        {
            return scalar(0);
        };

        virtual void correctBaffleLessFields()
        {}
        virtual void deformMesh()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
