/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2512                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2025 OpenCFD Ltd.          |
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
    approximation), expanded to include multi-group neutron diffusion,
    delayed neutron precursors and decay heat precursors for MSFR
    (Molten Salt Fast Reactor) applications, with complete coupling between
    flow, temperature and neutronics.

    Also includes Eulerian transport of dissolved fission product species
    (production, decay, wall deposition and associated decay heat flux), plus
    Lagrangian tracking of nuclear (fission product) particles undergoing
    decay, decay heat release and thermophoresis, and a surface film model
    for particle deposition.

    Converted from the native OpenFOAM application ``msfrParcelFoam``
    (see ``msfrParcelFoamNATIVE`` for the unmodified reference copy) into a
    run-time selectable foamForNuclear solver class.

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
#include "physicoChemicalConstants.H"
#include "mathematicalConstants.H"
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
    pimpleControl npimple_;

    //- Flow fields

    volVectorField U_;
    surfaceScalarField phi_;
    autoPtr<surfaceVectorField> Uf_;

    singlePhaseTransportModel laminarTransport_;

    //- Boussinesq / thermal / scalar-transport properties

    dimensionedScalar beta_;
    dimensionedScalar CpRef_;
    dimensionedScalar TRef_;
    dimensionedScalar Pr_;
    dimensionedScalar Prt_;
    dimensionedScalar Sc_;
    dimensionedScalar Sct_;

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

    scalar cumulativeContErr_;

    //- Neutronics / delayed-neutron / decay-heat / fission-product transport

    IOdictionary nuclearPropertiesDict_;

    //- Neutron transport constants
    label nEnergyGroups_;
    Switch useDoppler_;
    dimensionedScalar DopplerControl_;
    dimensionedScalar rhoRefXS_;
    dimensionedScalar TRefXS_;

    PtrList<dimensionedScalar> D_ref_;
    PtrList<dimensionedScalar> alphaD_ref_;
    PtrList<dimensionedScalar> Sa_ref_;
    PtrList<dimensionedScalar> alphaSa_ref_;
    PtrList<dimensionedScalar> Sf_ref_;
    PtrList<dimensionedScalar> alphaSf_ref_;
    PtrList<dimensionedScalar> Ss_ref_;
    PtrList<dimensionedScalar> alphaSs_ref_;
    PtrList<dimensionedScalar> Nu_;
    PtrList<dimensionedScalar> Ef_;
    PtrList<dimensionedScalar> invVel_;
    PtrList<dimensionedScalar> chiPrompt_;
    PtrList<dimensionedScalar> chiDelayed_;

    //- Delayed neutron precursors
    label nPrecGroups_;
    PtrList<dimensionedScalar> precLambda_;
    PtrList<dimensionedScalar> precBeta_;
    dimensionedScalar precBetaTot_;

    //- Decay heat precursors
    label nDecGroups_;
    PtrList<dimensionedScalar> decLambda_;
    PtrList<dimensionedScalar> decBeta_;
    dimensionedScalar decBetaTot_;

    //- Fission products (dissolved species transport, optional)
    bool fpTransport_;
    label nFPSpecies_;
    PtrList<dimensionedScalar> FPYield_;
    PtrList<dimensionedScalar> FPLambda_;
    PtrList<dimensionedScalar> FPDp_;
    PtrList<dimensionedScalar> FPRho_;
    PtrList<dimensionedScalar> FPPd_;

    //- Neutron flux
    PtrList<volScalarField> flux_;
    volScalarField fluxTot_;

    //- Precursors
    PtrList<volScalarField> prec_;
    PtrList<volScalarField> dec_;

    //- Fission products transport fields
    PtrList<volScalarField> Dfp_;
    PtrList<volScalarField> fp_;
    PtrList<volScalarField> fps_;
    PtrList<volScalarField> fpp_;
    PtrList<volScalarField> fpd_;
    PtrList<volScalarField> gamma_;
    PtrList<volScalarField> d_;
    PtrList<volScalarField> qdecay_wall_;

    //- Temperature/density-dependent nuclear properties
    PtrList<volScalarField> D_;
    PtrList<volScalarField> Sa_;
    PtrList<volScalarField> Sf_;
    PtrList<volScalarField> Ss_;
    PtrList<volScalarField> M_;

    autoPtr<volScalarField> diffCoeffAlbedo_;

    volScalarField unitField_;
    volScalarField zeroField_;

    volScalarField Rf_;
    volScalarField q_;
    volScalarField qPrompt_;
    volScalarField qDecay_;
    volScalarField Q_;
    volScalarField Qprompt_;
    volScalarField Qdecay_;

    volScalarField NuRf_;
    volScalarField EfRf_;

    PtrList<volScalarField> neutronSource_;

    volScalarField rhoNorm_;
    volScalarField logT_;

    //- Reactor control
    Switch tempCoupling_;
    Switch criticality_;
    dimensionedScalar KeffInit_;
    dimensionedScalar Qnominal_;

    volScalarField Keff_;
    volScalarField reactivity_;


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

        //- Solve thermal-hydraulics + neutronics
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
