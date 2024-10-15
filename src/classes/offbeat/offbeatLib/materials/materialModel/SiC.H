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
    Foam::SiC

Description
    SiC properties taken mostly from PARFUME;

    The material property and behavioral models can selected with the specific
    keywords. Otherwise, default models will be used.

    \htmlinclude materials_SiC.html

SourceFiles
    SiC.C

\mainauthor
    F. Xiang - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour) & XJTU (Xi'an JiaoTong
    University, China, Nuclear Thermal-hydraulic Laboratory)

\contribution
    A. Scolaro - EPFL
    E. Brunetto - EPFL

\date
    December 2022

\*---------------------------------------------------------------------------*/

#ifndef SiC_H
#define SiC_H

#include "materialModel.H"
#include "volFields.H"

#include "failureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class SiC Declaration
\*---------------------------------------------------------------------------*/

class SiC
:
    public materialModel
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        SiC(const SiC&);

        //- Disallow default bitwise assignment
        void operator=(const SiC&);


    //- Scaling parameters for UQ and SA

        //- Group name for parameter definition
        declareNameWithVirtual(SiC, group_, group);

        //- Density scaling factor [-]
        declareParameterWithVirtual(SiC, F_rho_, F_rho);

        //- Specific heat scaling factor [-]
        declareParameterWithVirtual(SiC, F_Cp_, F_Cp);

        //- Thermal conductivity scaling factor [-]
        declareParameterWithVirtual(SiC, F_k_, F_k);

        //- Emissivity scaling factor [-]
        declareParameterWithVirtual(SiC, F_emissivity_, F_emissivity);

        //- Elastic Young's modulus scaling factor [-]
        declareParameterWithVirtual(SiC, F_E_, F_E);

        //- Poisson's ratio scaling factor [-]
        declareParameterWithVirtual(SiC, F_nu_, F_nu);

        //- Themral expansion (alpha*DeltaT) scaling factor  [-]
        declareParameterWithVirtual(SiC, F_alphaT_, F_alphaT);

        // //- Density offset [kg/m3]
        declareParameterWithVirtual(SiC, delta_rho_, delta_rho);

        //- Specific heat offset [J/kg.K]
        declareParameterWithVirtual(SiC, delta_Cp_, delta_Cp);

        //- Thermal conductivity offset [W/m.K]
        declareParameterWithVirtual(SiC, delta_k_, delta_k);

        //- Emissivity offset [-]
        declareParameterWithVirtual(SiC, delta_emissivity_, delta_emissivity);

        //- Elastic Young's modulus offset [N/m2]
        declareParameterWithVirtual(SiC, delta_E_, delta_E);

        //- Poisson's ratio offset [-]
        declareParameterWithVirtual(SiC, delta_nu_, delta_nu);

        //- Themral expansion (alpha*DeltaT) offset [-]
        declareParameterWithVirtual(SiC, delta_alphaT_, delta_alphaT);


public:

    //- Runtime type information
    TypeName("SiC");


    // Constructors

        //- Construct from mesh, global options, materials, and dictionary
        SiC
        (
            const fvMesh& mesh,
            const dictionary& materialModelDict,
            const labelList& addr
        );


    //- Destructor
    ~SiC();


    // Access Functions


    // Member Functions

        //- Update the behavioral models
        virtual void correctBehavioralModels( const scalarField& T );

        //- Return name of the material
        virtual word name() const
        {
            return "SiC";
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
