/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "uniformDirectionalMultiMaterialInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "Random.H"
// #include "setWriter.H"
#ifdef OPENFOAMFOUNDATION
#include "helperFuns.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(uniformDirectionalMultiMaterialInterface, 0);
    addToRunTimeSelectionTable
    (
        multiMaterialInterface, 
        uniformDirectionalMultiMaterialInterface, 
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformDirectionalMultiMaterialInterface::uniformDirectionalMultiMaterialInterface
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    uniformMultiMaterialInterface(mesh, dict),
    rotate_(readBool(dict.lookup("rotated"))),
    rotation_(),
    defaultWeightEps_(dict.lookup("defaultWeights")),
    defaultWeightGrad_(defaultWeightEps_)
{
    if (rotate_)
    {
        rotation_ = coordinateRotation::New(dict.subDict("rotation"));
    }
    
    if (dict.found("defaultWeightsGrad"))
    {
        defaultWeightGrad_ = vector(dict.lookup("defaultWeightsGrad"));
    }
    
#ifdef OPENFOAMFOUNDATION
    if (debug)
    {
        autoPtr<MeshedSurface<face> > surf(helperFuns::meshedUnitSphere());
        surf->movePoints(surf->points() + vector(1,0,0));

        scalarField w(surf->size());
        const vectorField& nf = surf->faceNormals();
        
        if (rotate_)
        {
            rotation_->updatePoints(surf->faceCentres());
            w = mag(cmptMultiply(rotation_->invTransform(nf), defaultWeightEps_));
        }
        else
        {
            w = mag(cmptMultiply(nf, defaultWeightEps_));
        }

        helperFuns::writePatchFieldVTK
        (
            primitivePatch(SubList<face>(surf(), surf().size()), surf->points()), 
            w,
            "sphere",
            "weights",
            mesh.time().timeName()
        );
        Info<< "Directional weights written in VTK format as "
               """weights_sphere.vtk""" << nl << endl;
    }
#endif
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::uniformDirectionalMultiMaterialInterface::~uniformDirectionalMultiMaterialInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- Return a surfaceScalarField of interface weights
void Foam::uniformDirectionalMultiMaterialInterface::correct()
{
    Foam::uniformMultiMaterialInterface::correct();

    const fvMesh& referenceMesh = ( 
        mesh_.foundObject<fvMesh>("referenceMesh")
    ) ? (
        mesh_.lookupObject<fvMesh>("referenceMesh")
    ) : (
        mesh_
    );
        
    const vectorField& Sf = referenceMesh.Sf().field();
    const scalarField& magSf = referenceMesh.magSf().field();
    
    vectorField nf = Sf / magSf;
    
    if (rotate_)
    {        
    #ifdef OPENFOAMFOUNDATION
        rotation_->updatePoints(referenceMesh.Cf().internalField());
        nf = rotation_->invTransform(nf);
    #elif OPENFOAMESI
        // rotation_->updatePoints(surf->faceCentres());        
        nf = (rotation_->R().T() & nf);
    #endif
    }

    weightsEps_.field() = mag(cmptMultiply(nf,defaultWeightEps_));
    weightsGrad_.field() = mag(cmptMultiply(nf,defaultWeightGrad_));
}

// ************************************************************************* //
