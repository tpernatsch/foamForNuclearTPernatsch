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
    Foam::diffCoefModel

Description
    Mother class for fission products diffusion coefficient models.

SourceFiles
    diffCoefModel.C

\mainauthor
    F. Xiang - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour) & XJTU (Xi'an JiaoTong
    University, China, Nuclear Thermal-hydraulic Laboratory)

\contribution
    A. Scolaro - EPFL
    E. Brunetto - EPFL

\date
    June 2023

\*---------------------------------------------------------------------------*/

#ifndef diffCoefModel_H
#define diffCoefModel_H

#include "primitiveFields.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class diffCoefModel Declaration
\*---------------------------------------------------------------------------*/

class diffCoefModel
{
  // Private data


    // Private Member Functions

        //- Disallow default bitwise copy construct
        diffCoefModel(const diffCoefModel&);

        //- Disallow default bitwise assignment
        void operator=(const diffCoefModel&);

protected:



        //- Reference to object mesh
        const fvMesh& mesh_;

        // Protected data
        const bool crackFp_;
        const scalar failureTime_;

        bool isCracked() const
        {
            return crackFp_ || mesh_.time().value() >= failureTime_;
        }

        wordList FpNamesMat_;

        scalarList d1_;
        scalarList d2_;
        scalarList q1_;
        scalarList q2_;



public:

    //- Runtime type information
    TypeName("diffCoefModel");

    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        diffCoefModel,
        dictionary,
        (const fvMesh& mesh, const dictionary& dict, const word defaultModel),
        (mesh, dict, defaultModel)
    );

    // Constructors

        //- Construct from dictionary
        diffCoefModel
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word defaultModel
        );


    // Selectors

        //- Select from dictionary
        static autoPtr<diffCoefModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word defaultModel
        );

    //- Destructor
    virtual ~diffCoefModel();


    // Member Functions

    //- Update diffusion coefficient
    virtual void updateCoef(volScalarField& sf, const volScalarField& T, const labelList& addr, const word FpName) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
