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
    Foam::fuelMaterial

Description
    Parent class for all fuel materials (UO2, UPuO2 etc.).

    It defines some variables specific to fuels such as dishFraction_ or 
    enrichment_.

    It is not runTime selectable.

SourceFiles
    fuelMaterial.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef fuelMaterial_H
#define fuelMaterial_H

#include "materialModel.H"
#include "volFields.H"
#include "densificationModel.H"
#include "swellingModel.H"
#include "relocationModel.H"
#include "failureModel.H"
#include "poreVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fuelMaterial Declaration
\*---------------------------------------------------------------------------*/

class fuelMaterial
:
    public materialModel
{
    // Private data


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fuelMaterial(const fuelMaterial&);

        //- Disallow default bitwise assignment
        void operator=(const fuelMaterial&);

protected:

        //- densification model
        autoPtr<densificationModel> densification_;
        
        //- swelling model
        autoPtr<swellingModel> swelling_;
        
        //- relocation model
        autoPtr<relocationModel> relocation_;

        //- failure model
        autoPtr<failureModel> failure_;
        
        //- Pore velocity model
        autoPtr<poreVelocityModel> poreVelocity_;

        //- Initial enrichment
        const scalar enrichment_;

        //- Initial grain radius
        const scalar rGrain_;

        //- Initial oxygenMetal ratio
        const scalar oxygenMetalRatio_;

        //- Density Fraction
        const scalar densityFrac_;

        //- Porosity [-]
        volScalarField& porosity_;

        //- U composition - weight fractions
        //- namely: (U234, U235, U236, U238)
        const scalarList isotopesU_;

        //- Pu composition - weight fractions
        //- namely: (Pu238, Pu239, Pu240, Pu241, Pu242)
        const scalarList isotopesPu_;

        //- Am composition - weight fractions
        //- namely: (Am241, Am243)
        const scalarList isotopesAm_;

        // Ratio U/Metal in MOX fuel (%weight)
        const scalar ratioUMetal_;

        // Ratio Pu/Metal in MOX fuel (%weight)
        const scalar ratioPuMetal_;

        // Ratio Am/Metal in MOX fuel (%weight)
        const scalar ratioAmMetal_;

        //- Dish fraction
        scalar dishFraction_;

public:

    //- Runtime type information
    TypeName("fuelMaterial");


    // Constructors

        //- Construct from mesh, materials, and dictionary
        fuelMaterial
        (
            const fvMesh& mesh, 
            const dictionary& materialModelDict,
            const labelList& addr
        );


    //- Destructor
    ~fuelMaterial();
    

    // Access functions    

        //- Return the initial enrichment
        virtual scalar enrichment() const
        {
            return enrichment_;
        }

        //- Return the initial grain radius
        virtual scalar rGrain() const
        {
            return rGrain_;
        }

        //- Return the initial oxygen to metal ratio
        virtual scalar oxygenMetalRatio() const
        {
            return oxygenMetalRatio_;
        }

        //- Return the density fraction
        virtual scalar densityFrac() const
        {
            return densityFrac_;
        }

        //- Return the isotopesU
        virtual scalarList isotopesU() const
        {
            return isotopesU_;
        }

        //- Return the isotopesPu
        virtual scalarList isotopesPu() const
        {
            return isotopesPu_;
        }

        //- Return the isotopesAm
        virtual scalarList isotopesAm() const
        {
            return isotopesAm_;
        }

        //- Return the ratioUMetal
        virtual scalar ratioUMetal() const
        {
            return ratioUMetal_;
        }

        //- Return the ratioPuMetal
        virtual scalar ratioPuMetal() const
        {
            return ratioPuMetal_;
        }

        //- Return the ratioAmMetal
        virtual scalar ratioAmMetal() const
        {
            return ratioAmMetal_;
        }
        
        //- Return the material dishing fraction (0 by default)
        virtual const scalar& dishFraction() const
        {
            return dishFraction_;
        } 

        // //- Return pointer to pore velocity model
        // virtual autoPtr<poreVelocityModel> poreVelocity()
        // {
        //     return poreVelocity_;
        // }

    // Member Functions

        //- Return name of the material
        virtual word name() const
        {
            return "fuelMaterial";
        }; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
