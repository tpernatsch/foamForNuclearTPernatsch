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
    Foam::thermalExpansionMartinUPuO2

Description
    Class modelling thermal expansion of MOX fuel according to Martin 
    correlation.
    Source : "Thermophysical properties of MOX and UO2 Fuels including the 
              effects of irradiation" - S.G.Popov et al. - Oak Ridge National 
              Laboratory.

    https://info.ornl.gov/sites/publications/Files/Pub57523.pdf

SourceFiles
    thermalExpansionMartinUPuO2.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef thermalExpansionMartinUPuO2_H
#define thermalExpansionMartinUPuO2_H

#include "thermalExpansionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermalExpansionMartinUPuO2 Declaration
\*---------------------------------------------------------------------------*/

class thermalExpansionMartinUPuO2
:
    public thermalExpansionModel
{
    // Private data

        //- Parameters for the model
        scalar par1_;
        scalar par2_;
        scalar par3_;
        scalar par4_;
        scalar par5_;
        scalar par6_;
        scalar par7_;
        scalar par8_;

        //- Perturbation parameter
        scalar perturb;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermalExpansionMartinUPuO2(const thermalExpansionMartinUPuO2&);

        //- Disallow default bitwise assignment
        void operator=(const thermalExpansionMartinUPuO2&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UPuO2Martin");

    // Constructors

        //- Construct from dictionary
        thermalExpansionMartinUPuO2
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel   
        );

    //- Destructor
    virtual ~thermalExpansionMartinUPuO2();


    // Member Functions
    
        //- Update thermalExpansion  
        virtual void correct(symmTensorField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
