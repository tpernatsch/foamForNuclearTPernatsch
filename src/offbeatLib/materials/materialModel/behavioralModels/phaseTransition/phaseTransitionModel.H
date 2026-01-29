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
    Foam::phaseTransitionModel

Description
    Parent class for phase transition model. 

    The parent class typeName is "none". If selected in the solverDict the
    phase change is not considered, but the field 'betaFraction' is created and 
    added to registry.

    If the 'betaFraction' file is present in the starting time folder, 
    the internal field and boundary conditions are set according to the file. 

    Otherwise, the field is set to 0 in each cell, with zeroGradient BCs in all
    non-empty/-wedge patches.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        ...

        cladding
        {
            material Zircaloy;

            ...

            phaseTransitionModel     none;
        }

        ...
    }
    \endverbatim
    
SourceFiles
    phaseTransitionModel.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef phaseTransitionModel_H
#define phaseTransitionModel_H

#include "primitiveFields.H"
#include "fvMesh.H"
#include "volFields.H"

#include "userParameters.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class phaseTransitionModel Declaration
\*---------------------------------------------------------------------------*/

class phaseTransitionModel
{
    // Private data
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        phaseTransitionModel(const phaseTransitionModel&);

        //- Disallow default bitwise assignment
        void operator=(const phaseTransitionModel&);

protected:
    
    // Protected data
        
        //- Reference to object mesh
        const fvMesh& mesh_;

        //- Reference to beta phase fraction field
        volScalarField& betaFraction_;

        
    //- Scaling parameters for UQ and SA
        
        //- Group name for parameter definition
        declareNameWithVirtual(phaseTransitionModel, group_, group);
    
        //- Beta fration scaling factor [-]
        declareParameterWithVirtual(phaseTransitionModel, 
            F_betaFraction_, F_betaFraction);
        
        //- Beta fraction offset [-]
        declareParameterWithVirtual(phaseTransitionModel, 
            delta_betaFraction_, delta_betaFraction);

public:

    //- Runtime type information
    TypeName("none");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseTransitionModel,
            dictionary,
            (const fvMesh& mesh, const dictionary& dict, const dictionary& baseDict),
            (mesh, dict, baseDict)
        );


    // Constructors

        //- Construct from dictionary
        phaseTransitionModel
        (
            const fvMesh& mesh,
            const dictionary& dict ,
            const dictionary& baseDict   
        );


    // Selectors

        //- Select from dictionary
        static autoPtr<phaseTransitionModel> New
        (   
            const fvMesh& mesh,
            const dictionary& dict,
            const word defaultModel 
        );


    //- Destructor
    virtual ~phaseTransitionModel();


    //- Access functions

        //- Return beta fraction field
        virtual const volScalarField& betaFraction() const
        {
            return betaFraction_;
        };
    

    // Member Functions

        //- Update beta fraction 
        virtual void correct(const scalarField& T, const labelList& addr)
        {};    
                
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
