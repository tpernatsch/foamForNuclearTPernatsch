/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::layerAdditionRemovalPolyTopoChanger

Description
    Derived from polyTopoChanger:

        - The list of mesh modifiers defining the mesh dynamics is not read 
        from the "meshModifiers" IOobject. Instead, the mesh modifiers are 
        hard-coded to be of the layerAdditionRemoval type. 
        - For each mesh modifier, quantities such as min and max layer 
        thickness are provided on a mesh-dependent basis.

    TODO:
        - Allow user to set min/max thickness?
        - Allow other type of meshModifiers?

SourceFiles
    layerAdditionRemovalPolyTopoChanger.C

\*---------------------------------------------------------------------------*/

#ifndef layerAdditionRemovalPolyTopoChanger_H
#define layerAdditionRemovalPolyTopoChanger_H

#include "polyTopoChanger.H"
#include "layerAdditionRemoval.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
            Class layerAdditionRemovalPolyTopoChanger Declaration
\*---------------------------------------------------------------------------*/

class layerAdditionRemovalPolyTopoChanger
:
    public polyTopoChanger
{
    // Private Member Functions

    //- Switch to evaluate thickness layer from volume
    const Switch thicknessFromVolume_;

protected:

    // Protected data

public:

    //- Runtime type information
    TypeName("layerAdditionRemovalPolyTopoChanger");


    // Constructors

        //- Read constructor given IOobject and a polyMesh
        layerAdditionRemovalPolyTopoChanger(const IOobject&, polyMesh&);

        //- Read constructor for given polyMesh
        layerAdditionRemovalPolyTopoChanger
        (
            polyMesh& mesh, 
            const Switch& thicknessFromVolume = true
        );

        //- Disallow default bitwise copy construction
        layerAdditionRemovalPolyTopoChanger
        (
            const layerAdditionRemovalPolyTopoChanger&
        ) = delete;


    //- Destructor
    virtual ~layerAdditionRemovalPolyTopoChanger()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //