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
    Foam::corrosionModel

Description
    Base class for patch-specific corrosion models. 
    The parent class typeName is "fromLatestTime". If selected the corrosion and
    the oxide layer does not evolve.

SourceFiles
    corrosionModel.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef corrosionModel_H
#define corrosionModel_H

#include "primitiveFields.H"
#include "surfaceFields.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class corrosionModel Declaration
\*---------------------------------------------------------------------------*/

class corrosionModel
{
    // Private data
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        corrosionModel(const corrosionModel&);

        //- Disallow default bitwise assignment
        void operator=(const corrosionModel&);

protected:
    
    // Protected data
        
        //- Reference to object mesh
        const fvMesh& mesh_;
        
        //- Update mesh or not
        bool updateMesh_;

public:

    //- Runtime type information
    TypeName("fromLatestTime");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            corrosionModel,
            dictionary,
            (const fvMesh& mesh),
            (mesh)
        );


    // Constructors

        //- Construct from dictionary
        corrosionModel
        (
            const fvMesh& mesh
        );


    // Selectors

        //- Select from dictionary
        static autoPtr<corrosionModel> New
        (   
            const fvMesh& mesh,
            const word& type
        );


    // Destructors

        virtual ~corrosionModel();


    // Access functions
    
        //- Decide whether to update the mesh or not 
        virtual bool updateMesh()
        {
            return updateMesh_;
        };


    // Member Functions
    
        //- Update oxide thickness 
        virtual void correct
        (
            surfaceScalarField& oxideThickness,
            surfaceScalarField& DOxideThickness,
            const label& patchID
        ) {};
    
        //- Update change in metal thickness 
        virtual void updateDMetalThickness
        (
            const surfaceScalarField& DOxideThickness,
            surfaceScalarField& DMetalThickness,
            const label& patchID
        ) {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
