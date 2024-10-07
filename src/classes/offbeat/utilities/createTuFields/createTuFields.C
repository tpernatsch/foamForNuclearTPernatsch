/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

Application
    createTuFields

Description
    This application is intended to be used for the coupling with the 
    TRANSURANUS 1.5D fuel performance code. It simply reads TU fields from a 
    dictionary located in "$FOAM_CASE/system/TuOFFBEATdict" and creates OFFBEAT
    fields performing due interpolations on the OFFBEAT mesh.

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro - EPFL\n

\date 
    September 2023

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "HashTable.H"
#include "InterpolateTables.H"
#include "scalarFieldFieldINew.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Define TuField class
    #include "TuFieldClass.H"

    // Include hash tables with legend of Tu-to-OFFBEAT field names, units, etc
    #include "Tu2OFFBEATtable.H"

    IOdictionary TuOFFBEATdict
    (
        IOobject
        (
            "TuOFFBEATdict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    // Read time from dictionary
    scalar time(readScalar(TuOFFBEATdict.lookup("time")));

    // Set runTime to this value
    runTime.setTime(time, 0);

    // Get the names of the Tu fields in dict
    List<keyType> TuFieldNames(TuOFFBEATdict.subDict("fields").keys());

    // Create subDict of fields
    dictionary fieldsDict(TuOFFBEATdict.subDict("fields"));

    // Get Tu radial and axial points
    scalarField rValues(TuOFFBEATdict.lookup("radialPositions"));
    scalarField zValues(TuOFFBEATdict.lookup("axialPositions"));

    // Convert to meters
    rValues /= 1e3;
    zValues /= 1e3;

    // Initialize list of Tu fields with no corresponding OFFBEAT field
    wordList notTranslatedFields(0);

    // Interpolate fields from TU to OFFBEAT mesh
    #include "interpolateTuFields.H"

    // Build cylindrical fields from cylindrical components
    #include "buildCylindricalFields.H"

    // Build cartesian field from cylindrical field
    #include "buildCartesianFields.H"

    // Print warning for not translated fields, if any
    if ( notTranslatedFields.size() > 0 )
    {
        Info << nl << "WARNING! The following Tu fields have not OFFBEAT "
             << "corresponding fields: " << nl << endl;
        
        forAll(notTranslatedFields, i)
        {
            Info << "  '" << notTranslatedFields[i] << "' " << endl;
        }
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
