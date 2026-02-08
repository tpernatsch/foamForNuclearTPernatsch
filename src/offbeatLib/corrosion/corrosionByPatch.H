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
    Foam::corrosionByPatch

Description
    A corrosion model is created for each patch listed in the
    corrosionOptions dictionary.

    A polyTopoChanger is created with layerAdditionRemoval polyMeshModifiers
    for each patch is created.

    At each iteration:
        - the patch corrosion models are corrected
        - the patch face points are moved with the metal thickness change
        - if it is the first iteration of the current time step, the
          polyTopoChanger checks whether the topology should be updated by
          removing or adding a cell layer.

    TODO:
    -   Allow corrosion models on faceZones instead than patches?

Usage
    In solverDict file:
    \verbatim
    corrosion byPatch;

    corrosionOptions
    {
        // The updateMesh switch is set to true by default
        updateMesh     true;

        // The subdictionary must contain at least one patch
        patches
        {
            cladOuter
            {
                corrosionModel          zircaloyOuterCorrosion;
                oxidationKineticsModel  EPRI-KWU-CE|Cathcart-Pawel;
            }
        }
    }
    \endverbatim

SourceFiles
    corrosionByPatch.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date
    March 2022

\*---------------------------------------------------------------------------*/

#ifndef corrosionByPatch_H
#define corrosionByPatch_H

#include "corrosion.H"
#include "volFields.H"
#include "Table.H"
#include "layerAdditionRemovalPolyTopoChanger.H"
#include "polyTopoChange.H"
#include "corrosionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class corrosionByPatch Declaration
\*---------------------------------------------------------------------------*/

class corrosionByPatch
:
    public corrosion
{
    // Private data

        //- PolyTopoChanger
        layerAdditionRemovalPolyTopoChanger polyTopoChanger_;

        //- PtrList of corrosion models per patch ("fromLatestTime" as default)
        PtrList<corrosionModel> corrosionModels_;

        //- List of addressing for corrosion models
        labelList patchIDList_;

        //- Current time index
        int currentTimeIndex_;

        //- Current outer iteration index
        int nOuterIter_;

protected:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        corrosionByPatch(const corrosionByPatch&);

        //- Disallow default bitwise assignment
        void operator=(const corrosionByPatch&);

public:

    //- Runtime type information
    TypeName("byPatch");


    // Constructors

        //- Construct from mesh, materials and dict
        corrosionByPatch
        (
            fvMesh& mesh,
            const materials& mat,
            const dictionary& corrosionDict
        );


    //- Destructor
    ~corrosionByPatch();


    // Member Functions

        //- Update the corrosion
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
