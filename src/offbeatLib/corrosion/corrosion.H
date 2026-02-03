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
    Foam::corrosion

Description
    Base class for handling corrosion models. The user can select the corrosion 
    model specifying the keyword "corrosion" in the solver dictionary 
    (solverDict).

    The parent class typeName is "constant". If selected in the solverDict
    the corrosion does not evolve. However, the two corrosion fields 
    ('oxideThickness' and 'DOxideThickness') are created, added to the registry 
    and written at the end of the time step. 

    If the files are present in the starting time folder, the internal field 
    and boundary conditions are set according to the file.

    Otherwise, the fields are set to 0 m in each cell, with zeroGradient BCs 
    in all non-empty/-wedge patches.

    \warning
    For backward-compatibility, `fromLatestTime` can be still used as the typename
    for the `constant` corrosion model. However, its use is deprecated. **The class
    needs to be updated so that the `constant` class actually does not create oxide fields

Usage
    In solverDict file:
    \verbatim
    corrosion constant;

    corrosionOptions
    {
        // The updateMesh switch is set to true by default
        updateMesh     true;
    }
    \endverbatim

SourceFiles
    corrosion.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    March 2022

\*---------------------------------------------------------------------------*/

#ifndef corrosion_H
#define corrosion_H

#include "primitiveFields.H"
#include "volFields.H"
#include "materials.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class corrosion Declaration
\*---------------------------------------------------------------------------*/

class corrosion
{
    // Private data
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        corrosion(const corrosion&);

        //- Disallow default bitwise assignment
        void operator=(const corrosion&);

protected:

    // Protected data
        
        //- Reference to mesh
        fvMesh& mesh_;

        //- Reference to the materials class
        const materials& mat_;

        //- Corrosion options dictionary
        const dictionary corrosionDict_;

        //- Oxide layer thickness
        surfaceScalarField oxideThickness_;

        //- Increment in oxide layer thickness
        surfaceScalarField DOxideThickness_;

        //- Change in base material (metal) thickness
        //  (not necessarily identical to DOxideThickness_)
        surfaceScalarField DMetalThickness_;

        //- Activate update of mesh with oxide thickness
        Switch updateMesh_;


public:

    //- Runtime type information
    TypeName("constant");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            corrosion,
            dictionary,
            (
                fvMesh& mesh,
                const materials& mat,
                const dictionary& corrosionDict
            ),
            (mesh, mat, corrosionDict)
        );


    // Constructors

        //- Construct from mesh, materials and dictionary
        corrosion
        (
            fvMesh& mesh,
            const materials& mat,
            const dictionary& corrosionDict
        );


    // Selectors

        //- Select from dictionary
        static autoPtr<corrosion> New
        (
            fvMesh& mesh,
            const materials& mat,
            const dictionary& solverDict
        );


    //- Destructor
    virtual ~corrosion();


    // Access Functions

        //- Return reference to mesh
        const fvMesh& mesh()
        {
            return mesh_;
        };

        //- Return oxideThickness
        const surfaceScalarField& oxideThickness()
        {
            return oxideThickness_;
        };

        //- Return oxideThickness increment
        const surfaceScalarField& DOxideThickness()
        {
            return DOxideThickness_;
        };


    // Member Functions
    
        //- Update the corrosion   
        virtual void correct(){};
    
        //- Move mesh according to oxide layer thickness
        virtual void updateMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
