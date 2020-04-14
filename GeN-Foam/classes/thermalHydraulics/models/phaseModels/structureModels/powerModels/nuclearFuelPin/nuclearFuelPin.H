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
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)
    Carlo Fiorina (carlo.fiorina@outlook.com / carlo.fiorina@epfl.ch)
    
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

    //- Power density of the fuel
    volScalarField powerDensity_;

    //- Fields representing pin centerline, outer fuel, inner cladding and 
    //  outer cladding temperatures
    volScalarField Tfi_;
    volScalarField Tfo_;
    volScalarField Tci_;
    volScalarField Tco_;

    //- Average fuel and cladding temperatures
    volScalarField Tfav_;
    volScalarField Tcav_;

    //- Scalars that are input in the IOdictionary for passing min/max fuel
    //  clad temperatures
    scalar Tfmax_;
    scalar Tfmin_;
    scalar Tcmax_;
    scalar Tcmin_;

    //- Number of cells in the 1-D radial fuel mesh
    labelList fuelMeshSize_;

    //- Number of cells in the 1-D radial clad mesh
    labelList cladMeshSize_;

    //- Sum of fuelMeshSize and cladMeshSize
    labelList meshSize_;

    //- Radial position of the mesh nodes
    scalarListList r_;

    //- Fuel inner radius
    scalarList rfi_;

    //- Fuel outer radius
    scalarList rfo_;

    //- Cladding inner radius
    scalarList rci_;

    //- Cladding outer radius
    scalarList rco_;

    //- Delta radius in fuel (i.e. mesh cell length)
    scalarList drf_;

    //- Delta radius in cladding (i.e. mesh cell length)
    scalarList drc_;

    //- Gap width
    scalarList drg_;

    //- Fuel volumetric heat capacity
    scalarList rhoCpf_;

    //- Cladding volumetric heat capacity
    scalarList rhoCpc_;

    //- Fuel thermal conductivity
    scalarList kf_;

    //- Cladding thermal conductivity
    scalarList kc_;

    //- Gap heat transfer coefficient
    scalarList gapH_;

    //- True if fuel has hole
    boolList hollowFuel_;

    //- List that maps from cell index to region index
    labelList cellToRegion_;

public:

    //- Runtime type information
    TypeName("nuclearFuelPin");

    // Constructors

        nuclearFuelPin
        (
            const structureModel& structure,
            const dictionary& dicts
        );

        //- Destructor
        virtual ~nuclearFuelPin();

    // Member Functions

        void updateLocalAvgGlobalMinMaxT
        (
            const label& starti,
            const label& endi,
            const scalarList& r,
            const scalar& dr,
            const scalarField& Trad,
            scalar& Tav,
            scalar& Tmin,
            scalar& Tmax
        );

        //- Update temperature profiles Tc, Tf, in a single domain cell with
        //  current system conditions
        void updateLocalTemperatureProfile
        (
            const label& celli, 
            const scalar& ht, 
            const scalar& h
        );

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

        //- Return surface temperature (i.e. outer cladding temperature)
        virtual void correctT(volScalarField& T) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermoSubStructureModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //