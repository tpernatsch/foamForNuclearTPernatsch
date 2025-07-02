/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2506                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2025 OpenCFD Ltd.          |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::meshHandler

Description
    Class to hold all the region meshes and perform the volumetric mappings

Usage
    In multiRegionCouplingDict
    \verbatim
        mappings
        {
            <Name of region ONTO which I am mapping>
            {
                <Name of region FROM which I am mapping>
                {
                    sourceFields    (nameOfSourceField1 nameOfSourceField2);
                    targetFields    (nameOfTargetField1 nameOfTargetField2);
                }
            }
            neutroRegion
            {
                fluidRegion
                {
                    sourceFields    ( T U mu alpha );
                    targetFields    (TCool U mu alpha );
                }
            }
        }
    \endverbatim

SourceFiles
    meshHandler.C

\*---------------------------------------------------------------------------*/

#ifndef mapping_H
#define mapping_H

#include "solver.H"
#include "meshToMesh.H"
#include "zeroGradientFvPatchFields.H"
#include "meshHandler.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class meshHandler Declaration
\*---------------------------------------------------------------------------*/

class meshHandler
{
    // Private Member Data

        //- Region meshes

        PtrList<fvMesh> meshes_;

        PtrList<fvMesh> mappingMeshes_;

        dictionary mappingDict_;
        PtrList<PtrList<meshToMesh>> mappingList_;
        PtrList<PtrList<volScalarField>> scalarCouplingFields_;
        PtrList<PtrList<volVectorField>> vectorCouplingFields_;
        List<List<void*>> couplingFields_;


public:

    //Outer loop controls


    // Constructors

        //- Construct from components
        meshHandler(const Time& runTime);

        //- Disallow default bitwise copy construction
        meshHandler(const meshHandler&) = delete;


    //- Destructor
    ~meshHandler();


    // Member Functions

        void mapAllFields(const Time& runTime);

        void mapTheseFields(const Time& runTime, wordList meshToMap);

        void initializeMappedFields
        (
            const Time& runTime,
            const bool initializeMappedFields
        );

        void createBaffleLessMeshes(const Time& runTime);

        void createCouplingFields(const Time& runTime);

        fvMesh& returnMesh(word meshName);
        fvMesh& returnMappingMesh(word meshName);
        bool contains(wordList list, word thisWord);

        template<class Type>
        void map
        (
            word srcFieldName,
            word tgtFieldName,
            label tgtFieldLabel,
            label srcFieldLabel,
            bool removeBaffles
        );

        template<class Type>
        void mapAndWrite
        (
            word srcFieldName,
            word tgtFieldName,
            label tgtFieldLabel,
            label srcFieldLabel,
            bool removeBaffles,
            bool initializeMappedFields
        );

        template<class Type>
        void createBaffleLessField
        (
            word fieldName,
            label regionLabel,
            fvMesh& baffleLessMesh,
            label& nCouplingFields
        );


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "meshHandlerTemplates.C"
#endif

#endif

// ************************************************************************* //
