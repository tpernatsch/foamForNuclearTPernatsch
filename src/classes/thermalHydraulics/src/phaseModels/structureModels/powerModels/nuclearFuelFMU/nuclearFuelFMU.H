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
    Foam::powerModels::nuclearFuelFMU

Description
    Model for representing a nuclear fuel from an external FMU. The current
    model relies on the exchange of 3 axially sampled information:
    - Surface temperature of the cladding (GeN-Foam FMI output)
    - Power density in the fuel (GeN-Foam FMI output)
    - Heat flux between the structure and the fluid (GeN-Foam FMI input)

    The axial sampling can be shared with the FMUs.

    It is possible to sample between 1 and N points along a lines by
    providing x, y, and z coordinates. These samples are extracted in
    the mesh via the Tsurface.nuclearFuelFMU and powerDensityNeutronics fields.

    GeN-Foam receives from the FMU the surface power uniformly sampled at the
    boundary between the structure (e.g cladding) and the fluid. Then it is
    possible to reconstruct the Tsurface.nuclearFuelFMU field by interpolating
    these samples. The interpolation between the heat flux samples is
    performed using the Radial Basis Function interpolation algorithm
    (https://en.wikipedia.org/wiki/Radial_basis_function_interpolation). In case
    the model contains a single pin, it is possible to set a uniform XY heat
    flux distribution by setting the epsilon parameter to (0 0 nz).

    An interesting property of this powerModel is that the external fuel models
    can be in 1 or N FMUs, and coupled to 1 or M powerModels.

Usage
    \verbatim
    // Example for a 3 fuel pin model
    powerModel
    {
        type                        nuclearFuelFMU;
        volumetricArea              ${{2.0*$volFrac / $r}}; // 2.0 * volumeFrac / rCladOuter

        // Sampling position, in this case 3 sampling lines
        xPos                        (0.34   0.50   0.66);
        yPos                        (0.05   0.05   0.05);
        axialPowerInterpolationMethod   linear;
        axialLocationsNameToFMU     gfAxialLocations;
        axialLocations
        (
            0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4
        );

        // Radial Basis Function interpolation for heat flux reconstruction
        radialBasisFunctionMethod   gaussian; // gaussian or polyharmonicSpline or kriging
        epsilon                     (10 0 100); // Used only for Gaussian Radial Basis Functions interpolation

        // FMI ports
        avgPowerDensityNameToFMU
        (
            gfPowerDensity0 gfPowerDensity1 gfPowerDensity2
        );
        axialProfilePowerDensityNameToFMU
        (
            gfAxialProfPwrDens0 gfAxialProfPwrDens1 gfAxialProfPwrDens2
        );
        TstructNameToFMU
        (
            gfTstruct0 gfTstruct1 gfTstruct2
        );
        heatFluxNameFromFMU // or rhoCpdTdtNameFromFMU
        (
            gfSurfacePower0 gfSurfacePower1 gfSurfacePower2
        );
        initPowerDensity    100e6;
        initTstruct         600;
        initHeatFlux        1e5;
    }
    \endverbatim


    \verbatim
    // Example for a 7 fuel pin model in a lattice using the polyharmonic spline
    // interpolation method, fuel length is smaller than total sampling length
    powerModel
    {
        type                        nuclearFuelFMU;
        volumetricArea              ${{2.0*$volFrac / $r}}; // 2.0 * volumeFrac / rCladOuter

        // Sampling position, in this case 3 sampling lines
        xyPosLattice
        {
            origin          (0 0);
            type            hexagonal; // square; // hexagonal;
            nElements       (3 3);
            pitch           0.5;
            rotationAngle   0;
            lattice // 7 pins
            (
                e F F
                 F F F
                  F F e
            );
        }
        axialPowerInterpolationMethod   linear;
        axialLocationsNameToFMU     gfAxialLocations;
        axialLocations
        (
            0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4
        );
        // Optional fuel length incase the axial location is bigger than the
        // fuel, usually the case when simulating cladding longer than fuel
        fuelLength  0.8;

        // Radial Basis Function interpolation for heat flux reconstruction
        radialBasisFunctionMethod   polyharmonicSpline;

        // FMI ports
        avgPowerDensityNameToFMU
        (
            gfPowerDensity0 gfPowerDensity1 gfPowerDensity2 gfPowerDensity3
            gfPowerDensity4 gfPowerDensity5 gfPowerDensity6
        );
        axialProfilePowerDensityNameToFMU
        (
            gfAxialProfPwrDens0 gfAxialProfPwrDens1 gfAxialProfPwrDens2
            gfAxialProfPwrDens3 gfAxialProfPwrDens4 gfAxialProfPwrDens5
            gfAxialProfPwrDens6
        );
        TstructNameToFMU
        (
            gfTstruct0 gfTstruct1 gfTstruct2 gfTstruct3 gfTstruct4 gfTstruct5
            gfTstruct6
        );
        heatFluxNameFromFMU // or rhoCpdTdtNameFromFMU
        (
            gfSurfacePower0 gfSurfacePower1 gfSurfacePower2 gfSurfacePower3
            gfSurfacePower4 gfSurfacePower5 gfSurfacePower6
        );
        // Optional fuel and cladding axial samples
        TFuelNameFromFMU
        (
            gfTFuel0 gfTFuel1 gfTFuel2 gfTFuel3
            gfTFuel4 gfTFuel5 gfTFuel6
        );
        TCladNameFromFMU
        (
            gfTCladding0 gfTCladding1 gfTCladding2 gfTCladding3
            gfTCladding4 gfTCladding5 gfTCladding6
        );
        initPowerDensity    100e6;
        initTstruct         600;
        initHeatFlux        1e5;
    }
    \endverbatim

    where the entries mean:
    \table
        Property     | Description                        | Type | Req'd | Dflt
        type            | Type name: nuclearFuelFMU       | word |  yes  | -
        volumetricArea  | Volumetric area of the porous medium | scalar | yes | -
        xPos            | List of X positions of the sampling lines | scalarList | conditional | -
        yPos            | List of Y positions of the sampling lines | scalarList | conditional | -
        xyPosLattice    | XY lattice map dictionnary | dict | conditional | -
        axialPowerInterpolationMethod | Interpolation method for LHGR | word | no | linear
        axialLocationsNameToFMU | FMI port name of axial location list | word | yes | -
        axialLocations  | Axial position list | scalarList | yes | -
        radialBasisFunctionMethod | Radial Basis Function method name (gaussian, polyharmonicSpline, kriging) | word | yes | -
        epsilon         | Epsilon shaping factor for the Gaussian RBF | scalarList | conditional | -
        avgPowerDensityNameToFMU | List of FMI port names for LHGR output | wordList | yes | -
        axialProfilePowerDensityNameToFMU | List of FMI port names for power density normalized distribution output | wordList | yes | -
        TstructNameToFMU | List of FMI port names for Tstructure output | wordList | yes | -
        heatFluxNameFromFMU | List of FMI port names for heat flux input | wordList | conditional | -
        rhoCpdTdtNameFromFMU | List of FMI port names for rhoCpdTdt input | wordList | conditional | -
        TFuelNameFromFMU | List of FMI port names for fuel temperature input | wordList | no | -
        TCladNameFromFMU | List of FMI port names for cladding temperature input | wordList | no | -
        initPowerDensity | Initial power density | scalar | yes | -
        initTstruct     | Initial structure temperature | scalar | yes | -
        initHeatFlux    | Initial heat flux | scalar | conditional | -
        initRhoCpdTdt   | Initial rho*Cp*dT/dt | scalar | conditional | -
    \endtable

\mainauthor
    Thomas Guilbaud (thomas.guilbaud@epfl.ch)

SourceFiles
    nuclearFuelFMU.C

\*---------------------------------------------------------------------------*/

#ifndef nuclearFuelFMU_H
#define nuclearFuelFMU_H

#include "powerModel.H"
#include "IOFieldFields.H"
#include "InterpolateTablesGF.H"
#include "radialBasisFunctionInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace powerModels
{

/*---------------------------------------------------------------------------*\
                           Class nuclearFuelFMU Declaration
\*---------------------------------------------------------------------------*/

class nuclearFuelFMU
:
    public powerModel
{
protected:

    //- Dictionary from which FMI variables are written
    IOdictionary fmiState_;

    //- Average fuel temperature field for coupling
    volScalarField Tfav_;

    //- Average cladding temperature field for coupling
    volScalarField Tcav_;

    //- Fields representing the average surface temperature of the structure
    volScalarField Tsurface_;

    //- Heat flux FMI input flag
    List<bool> isHeatFluxInput_;

    //- Fraction of power from neutronics that goes to this powerModel
    //  This takes into account for instance cases where some of the
    //  heating goes to the moderator
    scalarList fractionOfPowerFromNeutronics_;

    //- If true, update avg TFuel and avg TClad fields for neutronics coupling
    List<bool> isFuelTemperatureFieldFromFMU_;

    //- X-position of sampling lines
    scalarListList xPos_;

    //- Y-position of sampling lines
    scalarListList yPos_;

    //- Axial location
    List<scalarField> axialLoc_;

    //- Optional fuel length used for average power density
    List<scalar> fuelLength_;

    //- Axial location interpolation method
    List<interpolateTableBaseGF::interpolationMethod> zMethod_;

    //- Radial basis function method name (polyharmonicSpline, gaussian)
    List<word> radialBasisFunctionMethod_;

    //- Shape parameter for the radial basis function. List of 3 values, one per
    //  direction. Only used by gaussian
    scalarListList epsilon_;

    //- Kriging variogram type (gaussian, spherical or exponential)
    List<Foam::radialBasisFunctionInterpolation::variogramType> krigingType_;

    //- Kriging parameter a, only used by the kriging option
    scalarList krigingA_;

    //- Kriging parameter c, only used by the kriging option
    scalarList krigingC_;

    //- Kriging parameter c0, only used by the kriging option
    scalarList krigingC0_;

    //- Kriging weights for cell, per region, only computed once, only used by
    //  the kriging option
    List<List<scalarList>> krigingWeight_;

    //- List that maps the region index to the region name
    wordList regionIndexToRegionName_;

    //- List that maps from cell index to region index
    labelList cellToRegion_;

    //- Inverse matrix of the Radial Basis Function operator for each cellZone
    List<SquareMatrix<scalar>> invRBFmatrix_;

public:

    //- Runtime type information
    TypeName("nuclearFuelFMU");

    // Constructors

        nuclearFuelFMU
        (
            structure& structure,
            const dictionary& dicts
        );

        //- Destructor
        virtual ~nuclearFuelFMU();


    // Member Functions

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

        //- Correct Tsurface field using the heat flux inputs from FMUs
        //  The heat flux is interpolated
        void correctHeatFluxInputsFromFMUs
        (
            const volScalarField& HT,
            const volScalarField& H,
            label regioni
        );

        //- Sample Tsurface and power density over lines and transfer data to
        //  the FMUs
        void correctInputsForFMUs(label regioni);

        //- Return surface temperature (i.e. outer cladding temperature)
        virtual void correctT(volScalarField& T) const;

private:

    // Private Member Functions

        //- Compute a recommended RBF epsilon parameter to avoid peaks in the
        //  RBF interpolation
        scalar computeEspilonForRBF(const scalarList list) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace powerModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //