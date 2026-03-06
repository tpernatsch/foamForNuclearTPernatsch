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
    Foam::Steel1515Ti

Description
    Class modelling the austenitic Steel 15-15 Ti cladding material.

SourceFiles
    Steel1515Ti.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef Steel1515Ti_H
#define Steel1515Ti_H

#include "materialModel.H"
#include "volFields.H"

#include "swellingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class Steel1515Ti Declaration
\*---------------------------------------------------------------------------*/

class Steel1515Ti
:
    public materialModel
{
    // Private data

        //- irradiationGrowth model
        autoPtr<swellingModel> swelling_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        Steel1515Ti(const Steel1515Ti&);

        //- Disallow default bitwise assignment
        void operator=(const Steel1515Ti&);


        //- Scaling parameters for UQ and SA
    
            //- Group name for parameter definition
            declareNameWithVirtual(Steel1515Ti, group_, group);
    
            //- Density scaling factor [-]
            declareParameterWithVirtual(Steel1515Ti, F_rho_, F_rho);
    
            //- Specific heat scaling factor [-]
            declareParameterWithVirtual(Steel1515Ti, F_Cp_, F_Cp);
    
            //- Thermal conductivity scaling factor [-]
            declareParameterWithVirtual(Steel1515Ti, F_k_, F_k);
    
            //- Emissivity scaling factor [-]
            declareParameterWithVirtual(Steel1515Ti, F_emissivity_, F_emissivity);
    
            //- Elastic Young's modulus scaling factor [-]
            declareParameterWithVirtual(Steel1515Ti, F_E_, F_E);
    
            //- Poisson's ratio scaling factor [-]
            declareParameterWithVirtual(Steel1515Ti, F_nu_, F_nu);
    
            //- Themral expansion (alpha*DeltaT) scaling factor  [-]
            declareParameterWithVirtual(Steel1515Ti, F_alphaT_, F_alphaT);
    
            // //- Density offset [kg/m3]
            declareParameterWithVirtual(Steel1515Ti, delta_rho_, delta_rho);
    
            //- Specific heat offset [J/kg.K]
            declareParameterWithVirtual(Steel1515Ti, delta_Cp_, delta_Cp);
    
            //- Thermal conductivity offset [W/m.K]
            declareParameterWithVirtual(Steel1515Ti, delta_k_, delta_k);
    
            //- Emissivity offset [-]
            declareParameterWithVirtual(Steel1515Ti, delta_emissivity_, delta_emissivity);
    
            //- Elastic Young's modulus offset [N/m2]
            declareParameterWithVirtual(Steel1515Ti, delta_E_, delta_E);
    
            //- Poisson's ratio offset [-]
            declareParameterWithVirtual(Steel1515Ti, delta_nu_, delta_nu);
    
            //- Themral expansion (alpha*DeltaT) offset [-]
            declareParameterWithVirtual(Steel1515Ti, delta_alphaT_, delta_alphaT);

public:

    //- Runtime type information
    TypeName("Steel1515Ti");


    // Constructors

        //- Construct from mesh, global options, materials, and dictionary
        Steel1515Ti
        (
            const fvMesh& mesh, 
            const dictionary& materialModelDict,
            const labelList& addr
        );


    //- Destructor
    ~Steel1515Ti();


    // Member Functions
    
        //- Update the behavioral models
        virtual void correctBehavioralModels( const scalarField& T );
        
        //- Return the tensor of additional strains
        virtual void additionalStrains( symmTensorField& sf ) const;

        //- Return name of the material
        virtual word name() const
        {
            return "Steel1515Ti";
        }; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
