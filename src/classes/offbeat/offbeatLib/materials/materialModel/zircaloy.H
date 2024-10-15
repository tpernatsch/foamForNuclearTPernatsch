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
    Foam::zircaloy

Description
    Zircaloy properties taken mostly from MATPRO, BISON manual and IAEA nuclear
    material handbook.

    The thermomechanical properties are taken from different sources.

    The material property and behavioral models can selected with the specific
    keywords. Otherwise, default models will be used.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        //- Name of cellZone with zircaloy material model
        cladding
        {
            material zircaloy;

            densityModel            ZyIAEA;
            heatCapacityModel       ZyIAEA;
            conductivityModel       ZyRELAP;
            emissivityModel         ZyConstant;
            YoungModulusModel       ZyMATPRO;
            PoissonRatioModel       ZyMATPRO;
            thermalExpansionModel   ZyMATPRO;

            swellingModel           growthBISONZy;
            phaseTransition         ZyStatic;
            failureModel            none;
        }
 
        ...

    }
    \endverbatim

SourceFiles
    zircaloy.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef zircaloy_H
#define zircaloy_H

#include "materialModel.H"
#include "volFields.H"

#include "swellingModel.H"
#include "phaseTransitionModel.H"

#include "failureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class zircaloy Declaration
\*---------------------------------------------------------------------------*/

class zircaloy
:
    public materialModel
{
    // Private data

        //- irradiationGrowth model
        autoPtr<swellingModel> swelling_;

        //- phase transition model
        autoPtr<phaseTransitionModel> phaseTransition_;

        //- failure model
        autoPtr<failureModel> failure_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        zircaloy(const zircaloy&);

        //- Disallow default bitwise assignment
        void operator=(const zircaloy&);

        
    //- Scaling parameters for UQ and SA
    
        //- Group name for parameter definition
        declareNameWithVirtual(zircaloy, group_, group);

        //- Density scaling factor [-]
        declareParameterWithVirtual(zircaloy, F_rho_, F_rho);
        
        //- Specific heat scaling factor [-]
        declareParameterWithVirtual(zircaloy, F_Cp_, F_Cp);
    
        //- Thermal conductivity scaling factor [-]
        declareParameterWithVirtual(zircaloy, F_k_, F_k);
    
        //- Emissivity scaling factor [-]
        declareParameterWithVirtual(zircaloy, F_emissivity_, F_emissivity);
    
        //- Elastic Young's modulus scaling factor [-]
        declareParameterWithVirtual(zircaloy, F_E_, F_E);
    
        //- Poisson's ratio scaling factor [-]
        declareParameterWithVirtual(zircaloy, F_nu_, F_nu);
        
        //- Themral expansion (alpha*DeltaT) scaling factor  [-]
        declareParameterWithVirtual(zircaloy, F_alphaT_, F_alphaT);

        // //- Density offset [kg/m3]
        declareParameterWithVirtual(zircaloy, delta_rho_, delta_rho);
        
        //- Specific heat offset [J/kg.K]
        declareParameterWithVirtual(zircaloy, delta_Cp_, delta_Cp);
    
        //- Thermal conductivity offset [W/m.K]
        declareParameterWithVirtual(zircaloy, delta_k_, delta_k);
    
        //- Emissivity offset [-]
        declareParameterWithVirtual(zircaloy, delta_emissivity_, delta_emissivity);
    
        //- Elastic Young's modulus offset [N/m2]
        declareParameterWithVirtual(zircaloy, delta_E_, delta_E);
    
        //- Poisson's ratio offset [-]
        declareParameterWithVirtual(zircaloy, delta_nu_, delta_nu);
        
        //- Themral expansion (alpha*DeltaT) offset [-]
        declareParameterWithVirtual(zircaloy, delta_alphaT_, delta_alphaT);


public:

    //- Runtime type information
    TypeName("zircaloy");


    // Constructors

        //- Construct from mesh, global options, materials, and dictionary
        zircaloy
        (
            const fvMesh& mesh, 
            const dictionary& materialModelDict,
            const labelList& addr
        );


    //- Destructor
    ~zircaloy();


    // Access Functions


    // Member Functions
    
        //- Update the behavioral models
        virtual void correctBehavioralModels
        (
            const scalarField& T
        );

        //- Check if failure occurred
        virtual void checkFailure();
        
        //- Return the tensor of additional strains
        virtual void additionalStrains( symmTensorField& sf ) const;  

        //- Return name of the material
        virtual word name() const
        {
            return "zircaloy";
        }; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
