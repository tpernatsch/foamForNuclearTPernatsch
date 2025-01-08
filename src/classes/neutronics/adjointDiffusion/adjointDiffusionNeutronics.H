/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2412                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2024 OpenCFD Ltd.          |
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
    Foam::adjointDiffusionNeutronics

Description
    Derived class for adjointDiffusion neutronics.

\mainauthor of this file (and associated .C or included .H files):
    Equations: Jun Shi <junshi@berkeley.edu>
    Coding: Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>, EPFL (Swizerland)

    The coding was not tested. Use at your own risk!

SourceFiles
    adjointDiffusionNeutronics.C

\*---------------------------------------------------------------------------*/

#ifndef adjointDiffusionNeutronics_H
#define adjointDiffusionNeutronics_H

#include "neutronics.H"
#include "XS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class adjointDiffusionNeutronics Declaration
\*---------------------------------------------------------------------------*/

class adjointDiffusionNeutronics
:
    public neutronics
{
protected:

    //- Cross sections
    XS xs_;

    //- Fields
    volScalarField Dalbedo_;
    PtrList<volScalarField> flux_;
    PtrList<volScalarField> fluxStar_;
    PtrList<volScalarField> prec_;//prec per total volume
    PtrList<volScalarField> precStar_;//prec per volume salt
    volScalarField fluxStarAlbedo_;
    volScalarField defaultFlux_;
    volScalarField defaultPrec_;

    //- Derived fields
    volScalarField scatteringSourceExtra_;

    //- Calc
    List<scalar> calculated_beta_eff;
    volScalarField adjointPrec_beta_;
    volScalarField adjointFlux_chiP_;
    volScalarField adjointFlux_chiD_;

private:

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        adjointDiffusionNeutronics(const adjointDiffusionNeutronics&);

        //- Disallow default bitwise assignment
        void operator=(const adjointDiffusionNeutronics&);


public:

    //- Runtime type information
        TypeName("adjointDiffusionNeutronics");


    // Constructors

        //- Construct from mesh, read data from IOdictionaries
        adjointDiffusionNeutronics(fvMesh& mesh);

    //- Destructor
        virtual ~adjointDiffusionNeutronics();


    // Member Functions

        //- Return fluxes
        virtual const PtrList<volScalarField> fluxes() const
        {
            return flux_;
        }

        //- Return precursors
        virtual const PtrList<volScalarField> precursors() const
        {
            return prec_;
        }

        //- Return power
        virtual scalar power() const
        {
            scalar pTot = 0;
            forAll(flux_, energyI)
            {
                pTot += fvc::domainIntegrate
                (
                    flux_[energyI] * xs_.sigmaPow()[energyI]
                ).value();
            }
            return pTot;
        }

        // virtual void getCouplingFieldRefs
        // (
        //     const objectRegistry& srcTH,
        //     const meshToMesh& neutroToFluid,
        //     const objectRegistry& srcTM,
        //     const meshToMesh& neutroToMech
        // );

        // virtual void interpolateCouplingFields
        // (
        //     const meshToMesh& neutroToFluid,
        //     const meshToMesh& neutroToMech
        // );

        //- Correct/update the properties
        virtual void correctPhysics();
        virtual void correctTightlyCoupledPhysics();

        virtual scalar maxDeltaT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} //End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
