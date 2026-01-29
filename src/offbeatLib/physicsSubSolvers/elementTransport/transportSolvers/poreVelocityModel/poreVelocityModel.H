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
    Foam::poreVelocityModel

Description
    Mother class for pore velocity models. If selected, sets the velocity to 
    zero.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        fuel
        {
            material UPuO2;

            ...

            poreVelocityModel     none;
            // Multiplicator for pore velocity
            // alphaPoreVelocity 1;
        }

        ...
    }
    \endverbatim

SourceFiles
    poreVelocityModel.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2022

\*---------------------------------------------------------------------------*/

#ifndef poreVelocityModel_H
#define poreVelocityModel_H

#include "primitiveFields.H"
#include "objectRegistry.H"
#include "volFields.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class poreVelocityModel Declaration
\*---------------------------------------------------------------------------*/

class poreVelocityModel
{
    // Private data
    
    // Private Member Functions

    //- Disallow default bitwise copy construct
    poreVelocityModel(const poreVelocityModel&);

    //- Disallow default bitwise assignment
    void operator=(const poreVelocityModel&);

protected:
    
    // Protected data
    
    //- Reference to mesh
    const fvMesh& mesh_;

    //- Input/options dictionary
    const dictionary materialModelDict_;

    //- Reference to velocity vector field
    volVectorField& poreVelocity_;

    //- Optional correction factor for pore velocity (-)
    const scalar alphaPoreVelocity_;

    Switch regularizeGradT_;
    scalar regularizationScale_;
    volVectorField regGradT_;

public:

    //- Runtime type information
    TypeName("none");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        poreVelocityModel,
        dictionary,
        (
            const fvMesh& mesh, 
            const dictionary& materialModelDict,
            const word defaultModel
        ),
        (mesh, materialModelDict, defaultModel)
    );

    // Constructors

        //- Construct from mesh, materials, and dictionary
        poreVelocityModel
        (
            const fvMesh& mesh, 
            const dictionary& materialModelDict,
            const word defaultModel
        );


    // Selectors

        //- Select from dictionary
        static autoPtr<poreVelocityModel> New
        (
            const fvMesh& mesh, 
            const dictionary& materialModelDict,
            const word defaultModel
        );

    // Destructor
    virtual ~poreVelocityModel();

    // Access functions

        //- Return reference to pore velocity field
        const volVectorField& poreVelocity() const 
        {
            return poreVelocity_;
        }

    // Member functions

        //- Update Velocity
        virtual void correct
        ( 
            const labelList& addr 
        ){};

        virtual const volVectorField regularizeGradT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
