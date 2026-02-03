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
    Foam::densificationModel

Description
    Parent class for densification model. 

    The parent class typeName is "none". If selected in the solverDict the
    densification is not considered, but the densification field 
    'epsilonDensification' is created and added to registry.

    If the 'epsilonDensification' file is present in the starting time folder, 
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

            densificationModel     none;
        }

        ...
    }
    \endverbatim

SourceFiles
    densificationModel.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef densificationModel_H
#define densificationModel_H

#include "primitiveFields.H"
#include "fvMesh.H"
#include "volFields.H"

#include "userParameters.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class densificationModel Declaration
\*---------------------------------------------------------------------------*/

class densificationModel
{
    // Private data

    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        densificationModel(const densificationModel&);

        //- Disallow default bitwise assignment
        void operator=(const densificationModel&);

protected:
    
    // Protected data
        
        //- Reference to object mesh
        const fvMesh& mesh_;

        //- Reference to densification strain field
        volScalarField& epsilonDensification_;

        
    // Scaling parameters for UQ and SA
        
        //- Group name for parameter definition
        declareNameWithVirtual(densificationModel, group_, group);
    
        //- Densification strain scaling factor [-]
        declareParameterWithVirtual(densificationModel, 
            F_epsilonDensification_, F_epsilonDensification);
        
        //- Densification strain offset [-]
        declareParameterWithVirtual(densificationModel, 
            delta_epsilonDensification_, delta_epsilonDensification);

public:

    //- Runtime type information
    TypeName("none");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            densificationModel,
            dictionary,
            (const fvMesh& mesh, const dictionary& dict, const dictionary& baseDict),
            (mesh, dict, baseDict)
        );


    // Constructors

        //- Construct from dictionary
        densificationModel
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict
        );


    // Selectors

        //- Select from dictionary
        static autoPtr<densificationModel> New
        (   
            const fvMesh& mesh,
            const dictionary& dict,
            const word defaultModel 
        );


    //- Destructor
    virtual ~densificationModel();


    //- Access functions

        //- Return densification strain field
        virtual const volScalarField& epsilonDensification() const
        {
            return epsilonDensification_;
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
