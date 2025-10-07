/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2212                                                  |
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
    Foam::SP3Neutronics

Description
    Subclass for SP3 neutronics [@FIORINA201524] [@Fiorina2015ApplicationCodes] [@FIORINA2017419].

<b>Reference publications</b>

<ol type="1">
    <li> Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk,
    GeN-Foam: a novel OpenFOAM® based multi-physics solver for 2D/3D transient
    analysis of nuclear reactors, Nuclear Engineering and Design, Volume 294,
    2015, Pages 24-37, ISSN 0029-5493, https://doi.org/10.1016/j.nucengdes.2015.05.035.

    <li> Carlo Fiorina, Konstantin Mikityuk, "Application of the new GeN-Foam multi-physics
    solver to the European Sodium Fast Reactor and verification against available codes",
    Proceedings of ICAPP 2015, May 03-06, 2015 - Nice (France), Paper 15226
</ol>

\mainauthor of this file (and associated .C or included .H files):
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>, EPFL (Switzerland)

SourceFiles
    SP3Neutronics.C

\*---------------------------------------------------------------------------*/

#ifndef SP3Neutronics_H
#define SP3Neutronics_H

#include "neutronics.H"
#include "XS.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class SP3Neutronics Declaration
\*---------------------------------------------------------------------------*/

class SP3Neutronics
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
    PtrList<volScalarField> fluxStar2_;
    PtrList<volScalarField> prec_;//prec per total volume
    PtrList<volScalarField> precStar_;//prec per volume salt
    volScalarField fluxStarAlbedo_;
    volScalarField defaultFlux_;
    volScalarField defaultFlux2_;
    volScalarField defaultPrec_;

    //- Derived fields
    volScalarField neutroSource_;
    volScalarField delayedNeutroSource_;
    volScalarField scatteringSourceExtra_;

private:

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        SP3Neutronics(const SP3Neutronics&);

        //- Disallow default bitwise assignment
        void operator=(const SP3Neutronics&);


public:

    //- Runtime type information
        TypeName("SP3Neutronics");


    // Constructors

        //- Construct from mesh, read data from IOdictionaries
        SP3Neutronics(dynamicFvMesh& mesh);

    //- Destructor
        virtual ~SP3Neutronics();


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

        //- Correct/update the properties
        virtual void correctPhysics();
        virtual void correctTightlyCoupledPhysics();

        virtual scalar maxDeltaT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
