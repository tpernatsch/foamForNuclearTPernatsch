/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::burnupLassmann

Description
    Class to handle the depletion of 14 selected nuclides through the numerical
    solution of simplified Bateman's equations. To use this module, a neutronics
    subsolver must be enabled to provide the neutron flux form factor.

    This module is based on the Lassmann algorithm "K. Lassmann, C. O'Carroll,
    J. van de Laar, C.T. Walker, The radial distribution of plutonium in high 
    burnup UO2 fuels,Journal of Nuclear Materials".

    \htmlinclude burnup_burnupLassmann.html

SourceFiles
    burnupLassmann.C

\mainauthor
    Y. Robert - INSA, Lyon\n
    A. Scolaro, E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE,
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef burnupLassmann_H
#define burnupLassmann_H

#include "volFields.H"
#include "objectRegistry.H"
#include "IOdictionary.H"

#include "fuelMaterial.H"
#include "sliceMapper.H"

#include "primitiveFields.H"
#include "fvMesh.H"
#include "surfaceFields.H"

#include "burnupFromPower.H"
#include "globalOptions.H"

#include "InterpolateTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class burnupLassmann Declaration
\*---------------------------------------------------------------------------*/

class burnupLassmann
:
    public burnupFromPower
{
    // Private data

        //- Power density form factor
        volScalarField formFactor_;

        //- Burnup in MWd/kgU
        volScalarField BukgU_;

        //- Number of elements for depletion
        const scalar Nelements_;

        //- Atomic densities for U235
        volScalarField N_U235_;

        //- Atomic densities for U236
        volScalarField N_U236_;

        //- Atomic densities for Np237
        volScalarField N_Np237_;

        //- Atomic densities for Pu238
        volScalarField N_Pu238_;

        //- Atomic densities for U238
        volScalarField N_U238_;

        //- Atomic densities for Pu239
        volScalarField N_Pu239_;

        //- Atomic densities for Pu240
        volScalarField N_Pu240_;

        //- Atomic densities for Pu241
        volScalarField N_Pu241_;

        //- Atomic densities for Pu242
        volScalarField N_Pu242_;

        //- Atomic densities for Am241
        volScalarField N_Am241_;

        //- Atomic densities for Am243
        volScalarField N_Am243_;

        //- Atomic densities for Cm242
        volScalarField N_Cm242_;

        //- Atomic densities for Cm243
        volScalarField N_Cm243_;

        //- Atomic densities for Cm244
        volScalarField N_Cm244_;

        //- Nuclides index for U235
        const int iU235_  ;

        //- Nuclides index for U236
        const int iU236_  ;

        //- Nuclides index for Np237
        const int iNp237_ ;

        //- Nuclides index for Pu238
        const int iPu238_ ;

        //- Nuclides index for U238
        const int iU238_  ;

        //- Nuclides index for Pu239
        const int iPu239_ ;

        //- Nuclides index for Pu240
        const int iPu240_ ;

        //- Nuclides index for Pu241
        const int iPu241_ ;

        //- Nuclides index for Pu242
        const int iPu242_ ;

        //- Nuclides index for Am241
        const int iAm241_ ;

        //- Nuclides index for Am243
        const int iAm243_ ;

        //- Nuclides index for Cm242
        const int iCm242_ ;

        //- Nuclides index for Cm243
        const int iCm243_ ;

        //- Nuclides index for Cm244
        const int iCm244_ ;

        //- Microscopic fission cross section
        scalarField sf_;

        //- Microscopic capture cross section
        scalarField sc_;

        //- Thermal absorption microscopic cross section
        scalarField sa_th_;

        //- XS for U238(n,2n)Np237 reaction
        scalar sn2n_;

        //- Pellet inner radius (for each cell)
        scalarField rIn_;

        //- Pellet outer radius (for each cell) 
        scalarField rOut_;

        //- U235 enrichement
        scalarField enrichment_;

        //- Fuel density fraction    
        scalarField densityFraction_;    

        //- Convergence criterion
        scalar precision_; 

        //- Mapper for averaging quantities in axial slices
        const sliceMapper* mapper_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        burnupLassmann(const burnupLassmann&);

        //- Disallow default bitwise assignment
        void operator=(const burnupLassmann&);

public:

    //- Runtime type information
    TypeName("Lassmann");

    // Constructors

    //- Construct from mesh, global Options, materials and dict
    burnupLassmann
    (
        const fvMesh& mesh,
        const materials& mat,
        const dictionary& burnupDict
    );

    //- Destructor
    ~burnupLassmann();

    // Access Functions

    // Member Functions
    
    //- Update the burnup
    virtual void correct();

    //- Return a const reference to SigmaA_ 
    virtual const scalarField SigmaA() const;

    //- Return a const reference to the diffusion coeff field 
    virtual const scalarField neutronD() const;

    //- Calculate inner and outer radii
    virtual void calcRinRout(const labelList& addr);

    //- Initialize XS and nuclear data
    virtual void initNuclearData(const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
