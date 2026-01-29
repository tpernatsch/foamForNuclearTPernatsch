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
    Foam::swellingModel

Description
    Parent class for swelling model. 

    The parent class typeName is "none". If selected in the solverDict the
    swelling is not considered, but the swelling field 'epsilonSwelling' is 
    created and added to registry.

    If the 'epsilonSwelling' file is present in the starting time folder, 
    the internal field and boundary conditions are set according to the file. 

    Otherwise, the field is set to 0 in each cell, with zeroGradient BCs in all
    non-empty/-wedge patches.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        fuel
        {
            material UO2;

            ...

            swellingModel     none;
        }

        ...
    }
    \endverbatim

SourceFiles
    swellingModel.C

\todo
    Axial swelling not necessarily along z

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef swellingModel_H
#define swellingModel_H

#include "primitiveFields.H"
#include "fvMesh.H"
#include "volFields.H"

#include "userParameters.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class swellingModel Declaration
\*---------------------------------------------------------------------------*/

class swellingModel
{
    // Private data
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        swellingModel(const swellingModel&);

        //- Disallow default bitwise assignment
        void operator=(const swellingModel&);

protected:
    
    // Protected data
        
        //- Reference to object mesh
        const fvMesh& mesh_;

        //- Pointer to swelling strain field
        volSymmTensorField& epsilonSwelling_;

        
    //- Scaling parameters for UQ and SA
        
        //- Group name for parameter definition
        declareNameWithVirtual(swellingModel, group_, group);
    
        //- Swelling strain scaling factor [-]
        declareParameterWithVirtual(swellingModel, 
            F_epsilonSwelling_, F_epsilonSwelling);
        
        //- Swelling strain offset [-]
        declareParameterWithVirtual(swellingModel, 
            delta_epsilonSwelling_, delta_epsilonSwelling);

public:

    //- Runtime type information
    TypeName("none");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            swellingModel,
            dictionary,
            (const fvMesh& mesh, const dictionary& dict, const dictionary& baseDict),
            (mesh, dict, baseDict)
        );


    // Constructors

        //- Construct from dictionary
        swellingModel
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict
        );


    // Selectors

        //- Select from dictionary
        static autoPtr<swellingModel> New
        (   
            const fvMesh& mesh,
            const dictionary& dict,
            const word defaultModel 
        );


    //- Destructor
    virtual ~swellingModel();


    //- Access functions

        //- Return swelling strain field
        virtual const volSymmTensorField& epsilonSwelling() const
        {
            return epsilonSwelling_;
        };


    // Member Functions
    
        //- Update heatCapacity 
        virtual void correct(const scalarField& T, const labelList& addr) 
        {};    
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
