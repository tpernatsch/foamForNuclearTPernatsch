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
    Foam::solvers::nuclearParcelFoam

Description
    Transient solver for buoyant, turbulent, incompressible flow (Boussinesq
    approximation), coupled to Lagrangian tracking of nuclear (fission
    product) particles undergoing decay, decay heat release and
    thermophoresis, plus a surface film model for particle deposition.

    Converted from the native OpenFOAM application ``nuclearParcelFoam``
    (see ``nuclearParcelFoamNATIVE`` for the unmodified reference copy)
    into a run-time selectable foamForNuclear solver class.

SourceFiles
    nuclearParcelFoam.C

\*---------------------------------------------------------------------------*/

#ifndef nuclearParcelFoam_H
#define nuclearParcelFoam_H

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
                      Class nuclearParcelFoam Declaration
\*---------------------------------------------------------------------------*/

class nuclearParcelFoam
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

    scalar cumulativeContErr_;


public:

    TypeName("nuclearParcelFoam");

    //- Constructors

        nuclearParcelFoam
        (
            dynamicFvMesh& mesh
        );

        //- Destructor
        virtual ~nuclearParcelFoam(){}

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
