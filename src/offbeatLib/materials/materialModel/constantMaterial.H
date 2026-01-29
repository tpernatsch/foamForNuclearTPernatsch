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
    Foam::constantMaterial

Description
    Class for the user manual input of the thermomechanical properties of the 
    simulated materials.

    All thermoMechanical property models are set to "constant".

SourceFiles
    constantMaterial.C

Usage
    In solverDict file:
    \verbatim
    materialProperties byZone;

    // List of materials, one per cellZone.
    materials
    {
        ...

        cladding
        {
            material constant;

            rho         rho     [1 -3 0 0 0]    6560;
            Cp          Cp      [0 2 -2 -1 0]   285;
            k           k       [1 1 3 -1 0]    21.5;
            emissivity emissivity [0 0 0 0 0]   0.808642;
            E           E       [1 -1 -2 0 0]   9.93e10;
            nu          nu      [0 0 0 0 0]     0.37;
            G           G       [1 -1 -2 0 0]   3.6241e+10;
            alpha       alpha   [0 0 0 0 0]     6e-6;
            Tref        Tref    [0 0 0 1 0]     300;

            ...
        }

        ....
    }
    \endverbatim

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef constantMaterial_H
#define constantMaterial_H

#include "materialModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class constantMaterial Declaration
\*---------------------------------------------------------------------------*/

class constantMaterial
:
    public materialModel
{
    // Private data
        

    // Private Member Functions

        //- Disallow default bitwise copy construct
        constantMaterial(const constantMaterial&);

        //- Disallow default bitwise assignment
        void operator=(const constantMaterial&);

public:

    //- Runtime type information
        TypeName("constant");


    // Constructors

        //- Construct from mesh, materials, and dictionary
        constantMaterial
        (
            const fvMesh& mesh, 
            const dictionary& materialModelDict,
            const labelList& addr
        );


    //- Destructor
    ~constantMaterial();


    // Member Functions     

        //- Return material name
        virtual word name() const
        {
            return "constant";
        }; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
