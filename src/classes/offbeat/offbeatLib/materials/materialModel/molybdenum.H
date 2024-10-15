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
    Foam::molybdenum

Description
    Properties for molybdenum, taken mostly from BISON manual.

    The material property and behavioral models can selected with the specific
    keywords. Otherwise, default models will be used.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        //- Name of cellZone with molybdenum material model
        molybdenumDisc
        {
            material molybdenum;

            densityModel            MoConstant;
            heatCapacityModel       Mo;
            conductivityModel       Mo;
            emissivityModel         MoConstant;
            YoungModulusModel       Mo;
            PoissonRatioModel       MoConstant;
            thermalExpansionModel   Mo;
            
            // Lacking better models we use FeCrAl model from BISON manual.
            swellingModel           FeCrAl;
        }
    }
    \endverbatim

SourceFiles
    molybdenum.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef molybdenum_H
#define molybdenum_H

#include "materialModel.H"
#include "volFields.H"

#include "swellingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class molybdenum Declaration
\*---------------------------------------------------------------------------*/

class molybdenum
:
    public materialModel
{
    // Private data
    
        //- Constant density [kg/m3]
        dimensionedScalar rho_;
        
        //- swelling model
        autoPtr<swellingModel> swelling_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        molybdenum(const molybdenum&);

        //- Disallow default bitwise assignment
        void operator=(const molybdenum&);

public:

    //- Runtime type information
    TypeName("molybdenum");


    // Constructors


    //- Construct from mesh, materials, and dictionary
    molybdenum
    (
        const fvMesh& mesh, 
        const dictionary& materialModelDict,
        const labelList& addr
    );


    //- Destructor
    ~molybdenum();


    // Member Functions
    
        //- Update the behavioral models
        virtual void correctBehavioralModels( const scalarField& T );

        //- Return the tensor of additional strains
        virtual void additionalStrains( symmTensorField& sf ) const;

        //- Return name of the material
        virtual word name() const
        {
            return "molybdenum";
        }; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
