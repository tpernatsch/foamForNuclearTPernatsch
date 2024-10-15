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
    Foam::relocationModel

Description
    Parent class for relocation model. 

    The parent class typeName is "none". If selected in the solverDict the
    relocation is not considered, but the relocation fields
    'epsilonRelocation' and 'epsilonRecoveredRelocation' are created and ad0ded 
    to registry.

    If the files are present in the starting time folder, the internal field and 
    boundary conditions are set according to the file. 

    Otherwise, the fields are set to 0 in each cell, with zeroGradient BCs in all
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

            relocationModel     none;
        }

        ...
    }
    \endverbatim
    
SourceFiles
    relocationModel.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef relocationModel_H
#define relocationModel_H

#include "primitiveFields.H"
#include "fvMesh.H"
#include "volFields.H"

#include "userParameters.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class relocationModel Declaration
\*---------------------------------------------------------------------------*/

class relocationModel
{
    // Private data
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        relocationModel(const relocationModel&);

        //- Disallow default bitwise assignment
        void operator=(const relocationModel&);

protected:
    
    // Protected data
        
        //- Reference to object mesh
        const fvMesh& mesh_;

        //- Reference to relocation strain field
        volScalarField& epsilonRelocation_;

        //- Reference to recovered relocation strain field
        volScalarField& epsilonRecoveredRelocation_;

        
    //- Scaling parameters for UQ and SA
        
        //- Group name for parameter definition
        declareNameWithVirtual(relocationModel, group_, group);
    
        //- Relocation strain scaling factor [-]
        declareParameterWithVirtual(relocationModel, 
            F_epsilonRelocation_, F_epsilonRelocation);

        //- Relocation recovery strain scaling factor [-]
        declareParameterWithVirtual(relocationModel, 
            F_epsilonRelocationRecovery_, F_epsilonRelocationRecovery);
        
        //- Relocation strain offset [-]
        declareParameterWithVirtual(relocationModel, 
            delta_epsilonRelocation_, delta_epsilonRelocation);
        
        //- Relocation recovery strain offset [-]
        declareParameterWithVirtual(relocationModel, 
            delta_epsilonRelocationRecovery_, delta_epsilonRelocationRecovery);

public:

    //- Runtime type information
    TypeName("none");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            relocationModel,
            dictionary,
            (const fvMesh& mesh, const dictionary& dict),
            (mesh, dict)
        );


    // Constructors

        //- Construct from dictionary
        relocationModel
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Selectors

        //- Select from dictionary
        static autoPtr<relocationModel> New
        (   
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel 
        );


    //- Destructor
    virtual ~relocationModel();


    //- Access functions

        //- Return relocation strain field
        virtual const volScalarField& epsilonRelocation() const
        {
            return epsilonRelocation_;
        };

        //- Return relocation recovery strain field
        virtual const volScalarField& epsilonRecoveredRelocation() const
        {
            return epsilonRecoveredRelocation_;
        };

    // Member Functions
    
        //- Update relocation strain 
        virtual void correct(const scalarField& T, const labelList& addr)
        {};    
        
        //- Apply user paramters to given field
        template<class T>
        void applyParameters
        (
            Field<T>& sf, 
            const labelList& addr, 
            const dimensionedScalar& F, 
            const dimensionedScalar& delta
        ) const
        {    
            forAll(addr, i)
            {
                label cellI(addr[i]);

                sf[cellI] = F.value()*sf[cellI] + delta.value()*pTraits<T>::one;
            }
        }; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
