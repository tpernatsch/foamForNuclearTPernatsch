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
    Foam::UO2

Description
    UO2 material class. 

    The thermomechanical properties are taken from different sources, as Bison 
    manual, Frapcon manual or Relap manual.

    The material property and behavioral models can selected with the specific
    keywords. Otherwise, default models will be used.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        //- Name of cellZone with UO2 material model
        fuel
        {
            material UO2;

            densityModel            UO2Constant;
            heatCapacityModel       UO2Matpro;
            conductivityModel       UO2Matpro;
            emissivityModel         UO2Relap;
            YoungModulusModel       UO2Matpro;
            PoissonRatioModel       UO2Constant;
            thermalExpansionModel   UO2Relap;

            densificationModel      UO2Frapcon;
            swellingModel           UO2Frapcon;
            relocationModel         UO2Frapcon;
            failureModel            none;
        }

        ...

    }
    \endverbatim

SourceFiles
    UO2.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef UO2_H
#define UO2_H

#include "fuelMaterial.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class UO2 Declaration
\*---------------------------------------------------------------------------*/

class UO2
:
    public fuelMaterial
{
    // Private data

        //- Initial enrichment (i.e. U235/(U235+238) .wt% )
        const scalar enrichment_;

    // Private Member Functions

        // Set initial enrichment
        const scalar setEnrichment();

        //- Disallow default bitwise copy construct
        UO2(const UO2&);

        //- Disallow default bitwise assignment
        void operator=(const UO2&);

        
    // Scaling parameters for UQ and SA
        
        //- Group name for parameter definition
        declareNameWithVirtual(UO2, group_, group);
    
        //- Density scaling factor [-]
        declareParameterWithVirtual(UO2, F_rho_, F_rho);
        
        //- Specific heat scaling factor [-]
        declareParameterWithVirtual(UO2, F_Cp_, F_Cp);
    
        //- Thermal conductivity scaling factor [-]
        declareParameterWithVirtual(UO2, F_k_, F_k);
    
        //- Emissivity scaling factor [-]
        declareParameterWithVirtual(UO2, F_emissivity_, F_emissivity);
    
        //- Elastic Young's modulus scaling factor [-]
        declareParameterWithVirtual(UO2, F_E_, F_E);
    
        //- Poisson's ratio scaling factor [-]
        declareParameterWithVirtual(UO2, F_nu_, F_nu);
        
        //- Themral expansion (alpha*DeltaT) scaling factor  [-]
        declareParameterWithVirtual(UO2, F_alphaT_, F_alphaT);

        // //- Density offset [kg/m3]
        declareParameterWithVirtual(UO2, delta_rho_, delta_rho);
        
        //- Specific heat offset [J/kg.K]
        declareParameterWithVirtual(UO2, delta_Cp_, delta_Cp);
    
        //- Thermal conductivity offset [W/m.K]
        declareParameterWithVirtual(UO2, delta_k_, delta_k);
    
        //- Emissivity offset [-]
        declareParameterWithVirtual(UO2, delta_emissivity_, delta_emissivity);
    
        //- Elastic Young's modulus offset [N/m2]
        declareParameterWithVirtual(UO2, delta_E_, delta_E);
    
        //- Poisson's ratio offset [-]
        declareParameterWithVirtual(UO2, delta_nu_, delta_nu);
        
        //- Themral expansion (alpha*DeltaT) offset [-]
        declareParameterWithVirtual(UO2, delta_alphaT_, delta_alphaT);

public:

    //- Runtime type information
    TypeName("UO2");


    // Constructors

        //- Construct from mesh, materials, and dictionary
        UO2
        (
            const fvMesh& mesh, 
            const dictionary& materialModelDict,
            const labelList& addr
        );


    //- Destructor
    ~UO2();


    // Access Functions

        //- Return the initial enrichment
        virtual const scalar enrichment() const
        {
            return enrichment_;
        }

    // Member Functions
    
        //- Update the behavioral models
        virtual void correctBehavioralModels( const scalarField& T );

        //- Check if failure occurred 
        virtual void checkFailure();
        
        //- Return the tensor of additional strains
        virtual void additionalStrains( symmTensorField& sf ) const;  

        //- Return name of the material
        virtual word name() const
        {
            return "UO2";
        }; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
