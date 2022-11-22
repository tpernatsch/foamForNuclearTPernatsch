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
    Foam::thermoSubstructures::lumpedNuclearStructure

Description
    Model for representing a lumped nuclear-heated structure characterized
    by a heat resistance from surface to average, and a heat resitance from
    average to max
    .

Authors
    Carlo Fiorina (carlo.fiorina@outlook.com / carlo.fiorina@epfl.ch)
    
SourceFiles
    lumpedNuclearStructure.C

\*---------------------------------------------------------------------------*/

#ifndef lumpedNuclearStructure_H
#define lumpedNuclearStructure_H

#include "powerModel.H"
#include "IOFieldFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace powerModels
{

/*---------------------------------------------------------------------------*\
                           Class lumpedNuclearStructure Declaration
\*---------------------------------------------------------------------------*/

class lumpedNuclearStructure
:
    public powerModel
{
protected:

    //- Field (over the global mesh) of scalarFields (over a 1-D mesh of size
    //  subMeshSize_) representing the temperatures in each node
    IOFieldField<Field, scalar> T_;

    //- Fields representing max and surface temp
    volScalarField Tmax_;
    volScalarField Tsurface_;

    //- Fraction of power from neutronics that goes to this powerModel
    //  This take sinto account for instance cases where some of the 
    //  heating goes to the moderator
    scalarList fractionOfPowerFromNeutronics_;

    //- Number of nodes
    labelList nodesNumber_;

    //- Node representing fuel for couling purposes
    labelList nodeFuel_;

    //- Node representing clad (or other structure) for coupling purposes
    labelList nodeClad_;

    //- Heat conductances between nodes (including first and last 
    //  "half" resitances)
    scalarListList Hs_;

    //- Volumetric heat capacity for each node
    scalarListList rhoCp_;

    //- Fraction of volume of the structure of each node
    scalarListList volFraction_;

    //- Fraction of total power from neutronics that goes to each node
    scalarListList qFraction_;

    //- List that maps from cell index to region index
    labelList cellToRegion_;

    //- List that maps the region index to the region name
    wordList regionIndexToRegionName_;


public:

    //- Runtime type information
    TypeName("lumpedNuclearStructure");

    // Constructors

        lumpedNuclearStructure
        (
            structure& structure,
            const dictionary& dicts
        );

        //- Destructor
        virtual ~lumpedNuclearStructure();

    // Member Functions

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

} // End namespace thermoSubstructures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //