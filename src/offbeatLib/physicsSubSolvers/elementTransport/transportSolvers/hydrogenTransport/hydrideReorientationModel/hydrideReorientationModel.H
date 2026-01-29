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
    Foam::hydrideReorientationModel

Description:
    Parent class for Hydride reorientation.
    The parent class typeName is "none". If selected in the solverDict the
    hydride reorientation is not considered.

Usage
    In solverDict:
    \verbatim
    elementTransportOptions
    {
        solvers ( HNGD );

        HSolverOptions
        {
            stressDrivenDiffusion   true;
            hydrideReorientationModel   none;

            materials
            {
                ...
            }
        }
    }
    \endverbatim

SourceFiles
    hydrideReorientationModel.C

\mainauthor
    L. Verma - Paul Scherrer Institut (PSI)

\contribution
    

\date 
    April 2024

\*---------------------------------------------------------------------------*/

#ifndef hydrideReorientationModel_H
#define hydrideReorientationModel_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fuelMaterial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class hydrideReorientationModel Declaration
\*---------------------------------------------------------------------------*/

class hydrideReorientationModel
{
    // Private data
            
    // Private Member Functions
       
        //- Disallow default bitwise copy construct
        hydrideReorientationModel(const hydrideReorientationModel&);

        //- Disallow default bitwise assignment
        void operator=(const hydrideReorientationModel&);
     
protected:
        
        //- Reference to mesh
        const fvMesh& mesh_;

        //- Input/options dictionary
        const dictionary lawDict_;

        //- Fraction of radial hydrides
        volScalarField& fracRHy_;

        //- Concentration of radial hydrides (wppm)
        volScalarField& radialHydrides_;


public:

    //- Runtime type information
    TypeName("hydrideReorientationModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            hydrideReorientationModel,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& lawDict
            ),
            (mesh, lawDict)
        );

    // Constructors

        //- Construct from mesh and dictionary
        hydrideReorientationModel
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );

    // Selectors

        //- Select from dictionary
        static autoPtr<hydrideReorientationModel> New
        (   
            const fvMesh& mesh,
            const dictionary& lawDict
        );
    
    //- Destructor
    virtual ~hydrideReorientationModel();

    // Access Functions

    // Member Functions
    
    //- Correct
    virtual void correctReorientation
    (
        const labelList& addr
    ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
