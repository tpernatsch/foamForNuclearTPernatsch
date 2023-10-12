/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2212                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2022 OpenCFD Ltd.         |
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
    - Surface power between the structure and the fluid (GeN-Foam FMI input)

    It is possible to uniformly sample between 1 and N points along a L lines by
    providing x, y, zmin and zmax coordinates. These samples are extracted in 
    the mesh via the Tsurface.nuclearFuelFMU and powerDensityNeutronics fields.

    GeN-Foam receives from the FMU the surface power uniformly sampled at the 
    boundary between the structure (e.g cladding) and the fluid. Then it is 
    possible to reconstruct the Tsurface.nuclearFuelFMU field by interpolating 
    these samples. The interpolation between the surface power samples is 
    performed using the Radial Basis Function interpolation algorithm 
    (https://en.wikipedia.org/wiki/Radial_basis_function_interpolation). In case 
    the model contains a single pin, it is possible to set a uniform XY surface 
    power distribution by setting the epsilon parameter to (0 0 nz).

    An interesting property of this powerModel is that the external fuel models 
    can be in 1 or N FMUs, and coupled to 1 or M powerModels.

Usage
    \verbatim
    // Example for a 3 fuel pin model
    powerModel
    {
        type            nuclearFuelFMU;
        volumetricArea  ${{2.0*$volFrac / $r}}; // 2.0 * volumeFrac / rCladOuter
        // Sampling position, in this case 3 sampling lines
        xPos            (0.34   0.50   0.66);
        yPos            (0.05   0.05   0.05);
        epsilon         (10 0 100); // Used for Radial Basis Functions interpolation
        // FMI ports
        avgPowerDensityNameToFMU
        (
            gfPowerDensity0 gfPowerDensity1 gfPowerDensity2
        );
        axialProfilePowerDensityNameToFMU
        (
            gfAxialProfilePowerDensity0 gfAxialProfilePowerDensity1 gfAxialProfilePowerDensity2
        );
        TstructNameToFMU
        (
            gfTstruct0 gfTstruct1 gfTstruct2
        );
        surfacePowerNameFromFMU
        (
            gfSurfacePower0 gfSurfacePower1 gfSurfacePower2
        );
        initPowerDensity        100e6;
        initTstruct             600;
        initSurfacePower        1e5;
    }
    \endverbatim

\mainauthor 
    Thomas Guilbaud (thomas.guilbaud@epfl.ch)
    
SourceFiles
    nuclearFuelFMU.C

\*---------------------------------------------------------------------------*/

#ifndef nuclearFuelFMU_H
#define nuclearFuelFMU_H

#include "powerModel.H"
#include "IOFieldFields.H"
#include "Function1.H"
#include "InterpolateTables.H"

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

    //- Average fuel and cladding temperatures
    // volScalarField Tfav_;
    // volScalarField Tcav_;

    //- Fields representing the average surface temperature of the structure
    volScalarField Tsurface_;

    //- Fraction of power from neutronics that goes to this powerModel
    //  This takes into account for instance cases where some of the 
    //  heating goes to the moderator
    scalarList fractionOfPowerFromNeutronics_;

    //- X-position of sampling lines
    scalarListList xPos_;

    //- Y-position of sampling lines
    scalarListList yPos_;

    //- Axial location
    List<scalarField> axialLoc_;

    //- Axial location interpolation method
    List<interpolateTableBase::interpolationMethod> zMethod_;

    //- Shape parameter for the radial basis function. List of 3 values, one per
    //  direction 
    scalarListList epsilon_;

    //- List that maps the region index to the region name
    wordList regionIndexToRegionName_;
    
    //- List that maps from cell index to region index
    labelList cellToRegion_;
    
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

        //- Return surface temperature (i.e. outer cladding temperature)
        virtual void correctT(volScalarField& T) const;


private:

    // Private Member Functions

        //- Weight calculation for the Radial Basis Function interpolation
        //  https://en.wikipedia.org/wiki/Radial_basis_function_interpolation
        scalarList solveRadialBasisFunction
        (
            const scalarList xList,
            const scalarList yList,
            const scalarList zList,
            const scalarList vList,
            const scalarList eps
        ) const;

        //- Compute the interpolation value based on the Radial Basis Function
        //  https://en.wikipedia.org/wiki/Radial_basis_function_interpolation
        scalar radialBasisFucntion
        (
            const scalarList w,
            const scalarList xList,
            const scalarList yList,
            const scalarList zList,
            const scalar x,
            const scalar y,
            const scalar z,
            const scalarList eps
        ) const;

        //- Stringify a list
        word stringify(const scalarList list) const;

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