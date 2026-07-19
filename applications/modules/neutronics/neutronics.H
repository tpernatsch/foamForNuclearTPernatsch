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
    Foam::neutronics

Group
    none

Description
    Base class for neutronic models.
    Reference publications [@FIORINA201524] [@FIORINA201515226].

Usage

    In the neutronicsProperties file, for keff calculation in a critical
    reactor:
    \verbatim
        eigenvalueNeutronics        true;
    \endverbatim

    For a subcritical reactor driven by an external neutron source (only for
    diffusionNeutronics sub-solver):
    \verbatim
        eigenvalueNeutronics        false;
        externalSourceNeutronics    true;
    \endverbatim


\mainauthor of this file (and associated .C or included .H files):
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>, EPFL (Swizerland)

SourceFiles
    neutronics.C


\*---------------------------------------------------------------------------*/

#ifndef neutronics_H
#define neutronics_H

#include "volFields.H"
#include "rhoThermo.H"
#include "fvMatrices.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "coordinateSystem.H"
#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "volPointInterpolation.H"
#include "meshToMesh.H"
#include "solver.H"
#include "timeProfile.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace solvers
{
/*---------------------------------------------------------------------------*\
                          Class neutronics Declaration
\*---------------------------------------------------------------------------*/

class neutronics
:
    public solver,
    public IOdictionary
{

protected:

    //- Dictionary from which keff, power are read
    IOdictionary reactorState_;

    //- k effective
    scalar keff_;

    //- Target power
    scalar pTarget_;

    //- Old total power
    scalar pTotOld_;

    //- Fuel power density
    volScalarField powerDensity_;

    //- Power density for for powerStructure, in case of liquid
    //  fuel, or for liquid, in case of solid fuel
    volScalarField secondaryPowerDensity_;

    //- One group averaged neutron flux
    volScalarField oneGroupFlux_;

    //- Ptrs to coupling fields (i.e. defined on their original mesh)
    // const volScalarField* TFuelOrig_;
    // const volScalarField* TCladOrig_;
    // const volScalarField* TCoolOrig_;
    // const volScalarField* rhoCoolOrig_;
    // const volScalarField* TStructOrig_;
    // const volScalarField* TStructMechOrig_;
    // const volVectorField* UOrig_;

    //- Total volume fraction of fluid(s)
    // const volScalarField* alphaOrig_;

    //- Thermal diffusivity in kg/m/s
    // const volScalarField* alphatOrig_;

    //- Molecular viscosity in kg/m/s (== Pa*s)
    // const volScalarField* muOrig_;

    //- Mapped fuel temperature field from thermal-hydraulics
    volScalarField TFuel_;

    //- Mapped cladding temperature field from thermal-hydraulics
    volScalarField TClad_;

    //- Mapped coolant temperature field from thermal-hydraulics
    volScalarField TCool_;

    //- Mapped coolant density field from thermal-hydraulics
    volScalarField rhoCool_;

    //- This is also used to compute TDriveline
    volScalarField TStruct_;

    //- Mapped structure temprature field from thermal-mechanics
    volScalarField TStructMech_;

    //- Also mapped, but these are autoPtr as these field are not needed at all
    //  unless doing liquidFuel. Thus, avoid wasting memory by only setting
    //  these fields if they are actually needed
    autoPtr<volVectorField> UPtr_;
    autoPtr<volScalarField> alphaPtr_;
    autoPtr<volScalarField> alphatPtr_;
    autoPtr<volScalarField> muPtr_;

    //- Fields derived from mapped fields
    autoPtr<surfaceScalarField> phiPtr_;
    autoPtr<volScalarField> diffCoeffPrecPtr_;

    //- Mesh displacement
    // volVectorField disp_;

    //- Initialize residual
    scalar initialResidual_;

    //- Solver residual
    scalar residual_;

    //- Eigenvalue flag
    bool eigenvalueNeutronics_;

    //- External source flag
    bool externalSourceNeutronics_;

    //- Liquid fuel (e.g. true for MSRs)
    const bool liquidFuel_;

    //-
    pointField originalPoints_;

    //- External source file
    IOdictionary externalSource_;

    //- Source modulation time profile
    timeProfile externalSourceModulationTimeProfile_;


private:

    //- Disallow default bitwise copy construct
    neutronics(const neutronics&);

    //- Disallow default bitwise assignment
    void operator=(const neutronics&);

public:

    //- Runtime type information
    TypeName("neutronics")

    //- Constructors

        //- Construct from mesh
        neutronics(dynamicFvMesh& mesh);

    //- Selectors

        //- Return a reference to the selected model
        // static autoPtr<neutronics> New
        // (
        //     fvMesh& mesh,
        //     bool multiRegion = false
        // );

    //- Destructor
        virtual ~neutronics(){}

    //- Member Functions

        //- Return fluxes
        virtual const PtrList<volScalarField> fluxes() const = 0;

        //- Return precursors
        virtual const PtrList<volScalarField> precursors() const = 0;

        //- Return power
        virtual scalar power() const = 0;

        virtual void correctBaffleLessFields();

        template<class Type>
        void correctBaffleLessField(word fieldName, fvMesh& baffleLessMesh);

        virtual void deformMesh();


        //- Required to remove "abstractedness" from IOdictionary
        virtual bool writeData(Ostream& os) const
        {
            return os.good();
        }


        //- Access

            virtual scalar keff() const
            {
                return keff_;
            }

            const volScalarField& powerDensity() const
            {
                return powerDensity_;
            }

            const volScalarField& secondaryPowerDensity() const
            {
                return secondaryPowerDensity_;
            }

            // const volVectorField& disp() const
            // {
            //     return disp_;
            // }

            bool eigenvalue() const
            {
                return eigenvalueNeutronics_;
            }

            scalar initialResidual() const
            {
                return initialResidual_;
            }

            virtual scalar getResidual()
            {
                return residual_;
            };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} //End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
