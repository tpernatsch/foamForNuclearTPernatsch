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

Class
    Foam::thermoSubStructureModels::nuclearFuelPin

Description
    Model for representing a nuclear fuel pin with constant material
    properties. Models fuel, gap, cladding and solves for the heat equation
    with a finite differences scheme. Only radial heat diffusion is considered.
    Derivation of the equations is provided at the end of the present file.

Authors
    Carlo Fiorina (carlo.fiorina@outlook.com / carlo.fiorina@epfl.ch)
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)
    
SourceFiles
    nuclearFuelPin.C

\*---------------------------------------------------------------------------*/

#ifndef nuclearFuelPin_H
#define nuclearFuelPin_H

#include "powerModel.H"
#include "IOFieldFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace powerModels
{

/*---------------------------------------------------------------------------*\
                           Class nuclearFuelPin Declaration
\*---------------------------------------------------------------------------*/

class nuclearFuelPin
:
    public powerModel
{
protected:

    //- Field (over the global mesh) of scalarFields (over a 1-D mesh of size
    //  subMeshSize_) representing the radial fuel temperature profile
    //  across fuel and cladding
    IOFieldField<Field, scalar> Trad_;

    //- Power density of the fuel (stored in fields_ table for I/O purposes)
    volScalarField powerDensity_;

    //- Fields representing pin centerline, outer fuel, inner cladding and 
    //  outer cladding temperatures (stored in fields_ table for I/O purposes)
    volScalarField Tfi_;
    volScalarField Tfo_;
    volScalarField Tci_;
    volScalarField Tco_;

    //- Number of cells in the 1-D radial fuel mesh
    label fuelSubMeshSize_;

    //- Number of cells in the 1-D radial clad mesh
    label cladSubMeshSize_;

    //- Total number of cells;
    label subMeshSize_;

    //- Radial position of the mesh nodes
    scalarList r_;

    //- Lenght of each 1-D fuel mesh cell
    scalar dFuelSubMesh_;

    //- Lenght of each 1-D clad mesh cell
    scalar dCladSubMesh_;

    //- Fuel inner radius
    scalar rfi_;

    //- Fuel outer radius
    scalar rfo_;

    //- Cladding inner radius
    scalar rci_;

    //- Cladding outer radius
    scalar rco_;

    //- Delta radius in fuel
    scalar drf_;

    //- Delta radius in cladding
    scalar drc_;

    //- Gap width
    scalar drg_;

    //- Fuel density
    scalar rhof_;

    //- Cladding density
    scalar rhoc_;

    //- Fuel thermal conductivity
    scalar kf_;

    //- Cladding thermal conductivity
    scalar kc_;

    //- Fuel heat capacity
    scalar Cpf_;

    //- Cladding heat capacity
    scalar Cpc_;

    //- Gap heat transfer coefficient
    scalar gapH_;

    //- True if fuel has hole
    const bool hollowFuel_;

public:

    //- Runtime type information
    TypeName("nuclearFuelPin");

    // Constructors

        nuclearFuelPin
        (
            const structureModel& structure,
            const dictionary& dict,
            const word name
        );

        //- Destructor
        virtual ~nuclearFuelPin();

    // Member Functions

        //- Update temperature profiles Tc, Tf, in a single domain cell with
        //  current system conditions
        void updateLocalTemperatureProfile(label celli, scalar ht, scalar h);

        //- Return surface temperature (i.e. outer cladding temperature)
        virtual tmp<volScalarField> T() const;

        //- In two-phase solvers :
        //  HT = H1*frac1*Tfluid1 + H2*frac2*Tfluid2
        //  H =  H1*frac1 + H2*frac2
        //- In a mono-phase solver :
        //  HT = H*Tfluid
        //  H = H
        virtual void correct
        (
            const volScalarField& HT,
            const volScalarField& H
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermoSubStructureModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //