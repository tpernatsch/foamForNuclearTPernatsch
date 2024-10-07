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
    Foam::failureModel

Description
    Mother class for materials failure criteria.

SourceFiles
    failureModel.C

\todo
    what if two+ materials have failureModel activated? Check what happens to 
    failure field.

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef failureModel_H
#define failureModel_H

#include "primitiveFields.H"
#include "fvMesh.H"
#include "volFields.H"
#include "sliceMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class failureModel Declaration
\*---------------------------------------------------------------------------*/

class failureModel
{
    // Private data
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        failureModel(const failureModel&);

        //- Disallow default bitwise assignment
        void operator=(const failureModel&);

protected:
    
    // Protected data

        //- Reference to object mesh
        const fvMesh& mesh_;

        //- Material dictionary
        const dictionary matDict_;

        //- Failure model dictionary
        const dictionary failureModelDict_;

        //- Pointer to field for burst flag (1 if burst occurred)
        //  The field is initialized only if failure actually occurs
        autoPtr<volScalarField> failedMaterial_;

        //- Switch for failure occurrence
        Switch failed_;

        //- Switch to break simulation if failure of a material occurs
        Switch stopIfFailed_;

        // - Mapper for averaging quantities in axial slices
        const sliceMapper* mapper_;

public:

    //- Runtime type information
    TypeName("none");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            failureModel,
            dictionary,
            (
                const fvMesh& mesh, 
                const dictionary& dict
            ),
            (mesh, dict)
        );

    // Constructors

        //- Construct from dictionary
        failureModel
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Selectors

        //- Select from dictionary
        static autoPtr<failureModel> New
        (   
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel 
        );


    //- Destructor
    virtual ~failureModel();


    //- Access functions


    // Member Functions

        //- Function called by daughter classes to initialize failedMaterial   
        virtual void initializeFailedMaterialField();

        //- Check failure occurrence 
        virtual void checkFailure(const labelList& addr);  

        //- Return true if material is failed (overwritten by daughter classes)
        virtual bool isFailed(const labelList& addr)
        {
            return false;
        };

        //- Return failure criterion name
        virtual word criterionName()
        {
            return typeName_();
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
