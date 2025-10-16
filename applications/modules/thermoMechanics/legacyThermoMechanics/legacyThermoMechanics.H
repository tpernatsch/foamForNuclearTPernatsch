/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2506                                                  |
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
    Foam::legacyThermoMechanics

Description
    Derived class for linear elastic thermoMechanics.

    Reference publications [@FIORINA201524] [@FIORINA201515226].

    Author of this file (and associated .C or included .H files):
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>, EPFL (Switzerland)

SourceFiles
    legacyThermoMechanics.C

\*---------------------------------------------------------------------------*/

#ifndef legacyThermoMechanics_H
#define legacyThermoMechanics_H

#include "thermoMechanics.H"
// #include "XS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class legacyThermoMechanics Declaration
\*---------------------------------------------------------------------------*/

class legacyThermoMechanics
:
    public thermoMechanics//class inheritance
{
protected:

    // Protected data

        //- Properties
        Switch solveDisplacement_;
        Switch solveThermo_;
        Switch planeStress_;
        Switch linkedFuel_;
        dictionary TMEntries_;
        label TMZoneNumber_;
        volScalarField rho_;
        volScalarField rhoE_;
        volScalarField nu_;
        volScalarField threeKalpha_;
        volScalarField DT_;
        volScalarField C_;
        volScalarField rhoK_;
        volScalarField alpha_;
        volScalarField E_;
        volScalarField mu_;
        volScalarField lambda_;
        volScalarField threeK_;


        //- Fields
        volVectorField disp_;
        volScalarField TStruct_;
        volScalarField TStructRef_;
        volSymmTensorField sigmaD_;
        volVectorField divSigmaExp_;
        volScalarField gapWidth;

        //- Mapped field
        volScalarField powerDensityNeutronics_;

        //- Derived from mapped field
        // Used to check where the TH is mapping TStructFromTH_
        volScalarField isPorous_;

        //- Control
        int nCorr_;
        scalar convergenceTolerance_;
        Switch compactNormalStress_;
        scalar residual_;

private:

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        legacyThermoMechanics(const legacyThermoMechanics&);

        //- Disallow default bitwise assignment
        void operator=(const legacyThermoMechanics&);


public:

    //- Runtime type information
        TypeName("legacyThermoMechanics");


    // Constructors

        //- Construct from mesh, read data from IOdictionaries
        legacyThermoMechanics(dynamicFvMesh& mesh);

    //- Destructor
        virtual ~legacyThermoMechanics();


    // Member Functions

        //- Substitutes the old getFields. Set local class references to global
        //  coupling fields, to be called only once after construction. The
        //  fields are set via lookup in the src objectRegistry (i.e.
        //  fluidMeshNB for the current state of GeN-Foam)
        // virtual void getCouplingFieldRefs
        // (
        //     const objectRegistry& srcTH,
        //     const objectRegistry& srcN
        // );

        // //- The old getFields, yet operates on the local refs without having
        // //  the need for all the fields to be passed. This effectively bypasses
        // //  the old limitation of getFields, i.e. if someone wanted to
        // //  implement a new thermoMechanical model that requires more fields
        // //  for the coupling, then they'd need to modify the getFields function
        // //  in ALL of the other thermoMechanical models too
        // virtual void interpolateCouplingFields
        // (
        //     const meshToMesh& mechToFluid,
        //     const meshToMesh& mechToNeutro
        // );

        //- Correct/update the properties
        virtual void correctPhysics();
        virtual void correctTightlyCoupledPhysics();
        virtual scalar maxDeltaT();
        virtual scalar getResidual()
        {
            return residual_;
        };

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
