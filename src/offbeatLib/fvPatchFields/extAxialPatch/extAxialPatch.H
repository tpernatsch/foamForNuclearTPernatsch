/*---------------------------------------------------------------------------*\
            Copyright (c) 2021, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the ECI4FOAM source code library, which is an
	unofficial extension to OpenFOAM.
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
    Foam::externalIOObject::extAxialPatch

Group
    grpFieldFunctionObjects

Description
    FMI function object to extract sampled list of values along a patch on the
    Z-axis.

Usage
    \verbatim
    surfacePowerFMI
    {
        type                extAxialPatch;
        fieldName           gradT;
        fieldType           vector;
        patchName           cladOuter;
        outputName          surfacePowerAxial;
        axialLocationsNameFromFMU   offbeatAxialLocations;
        initialData         21{0};
        initialAxialLocations
        (
            0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1
        );
        scale               21.5; // W/m/K, Thermal conductivity of Zircaloy
        relaxationFactor    0.2;
    }
    // or
    surfacePowerFMI
    {
        type                extAxialPatch;
        fieldName           gradT;
        patchName           cladOuter;
        outputName          surfacePowerAxial;
        axialLocationsNameFromFMU   offbeatAxialLocations;
        initialData         0;
        scale               21.5; // W/m/K, Thermal conductivity of Zircaloy
        relaxationFactor    0.2;
    }
    \endvarbatim

SourceFiles
    extAxialPatch.C

\mainauthor
    Thomas Guilbaud, EPFL/Transmutex SA, Switzerland

\date
    September 2023

\*---------------------------------------------------------------------------*/

#if defined __has_include
#  if __has_include(<commDataLayer.H>)
#    include <commDataLayer.H>
#    define isCommDataLayerIncluded
#  endif
#endif

#ifdef isCommDataLayerIncluded

#ifndef externalIOObject_extAxialPatch_H
#define externalIOObject_extAxialPatch_H

#include "externalIOObject.H"
#include "meshSearch.H"
#include "vector.H"
#include "fvMesh.H"
#include "globalOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace externalIOObject
{

/*---------------------------------------------------------------------------*\
                          Class extAxialPatch Declaration
\*---------------------------------------------------------------------------*/

class extAxialPatch
:
    public externalIOObject
{
    // Private Data

    //- FMI port name of the output
    word outputName_;

    //- FMI port name of the axial locations
    word axialLocationsNameFromFMU_;

    //- Field name to extract
    word fieldName_;

    //- Field type to extract
    word fieldType_;

    //- Patch name to extract along a line
    word patchName_;

    //- List of values samples along a line
    scalarList profileData_;

    //- Scaling factor
    scalar scale_;

    //- Relaxation factor
    scalar relax_;

    //- Mesh
    const fvMesh& mesh_;

    //- Pin axial direction
    vector pinDirection_;

    //- Dictionary from which FMI variables are written
    IOdictionary fmiState_;

public:

    //- Runtime type information
    TypeName("extAxialPatch");


    // Constructors

        //- Construct from Time and dictionary
        extAxialPatch
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        extAxialPatch(const extAxialPatch&) = delete;

        //- No copy assignment
        void operator=(const extAxialPatch&) = delete;


    //- Destructor
    virtual ~extAxialPatch() = default;


    // Member Functions

        //- Read the randomise data
        virtual bool read(const dictionary&);

        virtual bool execute();

        virtual bool executeScalarField();

        virtual bool executeVectorField();

        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace externalIOObject
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

#endif

// ************************************************************************* //
