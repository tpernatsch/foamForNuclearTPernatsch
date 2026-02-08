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

\*---------------------------------------------------------------------------*/

#if defined __has_include
#  if __has_include(<commDataLayer.H>)
#    include <commDataLayer.H>
#    define isCommDataLayerIncluded
#  endif
#endif

#ifdef isCommDataLayerIncluded

#include "addToRunTimeSelectionTable.H"
#include "extAxialPatch.H"
#include "commDataLayer.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace externalIOObject
{
    defineTypeNameAndDebug(extAxialPatch, 0);
    addToRunTimeSelectionTable(externalIOObject, extAxialPatch, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalIOObject::extAxialPatch::extAxialPatch
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    externalIOObject(name, runTime, dict),
    outputName_(),
    fieldName_(dict.get<word>("fieldName")),
    fieldType_(dict.get<word>("fieldType")),
    patchName_(dict.get<word>("patchName")),
    profileData_(),
    scale_(dict.lookupOrDefault<scalar>("scale", 1.0)),
    relax_(dict.lookupOrDefault<scalar>("relaxationFactor", 1.0)),
    mesh_
    (
        refCast<const fvMesh>
        (
            time_.lookupObject<objectRegistry>
            (
                dict.getOrDefault("region", polyMesh::defaultRegion)
            )
        )
    ),
    pinDirection_(),
    fmiState_
    (
        IOobject
        (
            name,
            mesh_.time().timeName(),
            "uniform/fmiState",
            mesh_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    )
{
    // Pin direction
    const globalOptions& globalOpt
    (
        mesh_.lookupObject<globalOptions>("globalOptions")
    );
    pinDirection_ = vector(globalOpt.pinDirection());

    // Initialize and first execution
    read(dict);
    execute();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::externalIOObject::extAxialPatch::read(const dictionary& dict)
{
    // Extract FMI port name
    if (!dict_.found("outputName"))
    {
        FatalError<< "outputName not found in dict"
            << abort(FatalError);
    }
    outputName_ = dict_.get<word>("outputName");

    if (!dict_.found("axialLocationsNameFromFMU"))
    {
        FatalError << "axialLocationsNameFromFMU not found in dict" << abort(FatalError);
    }
    axialLocationsNameFromFMU_ = dict_.get<word>("axialLocationsNameFromFMU");

    // Extract patch label
    label patchI = mesh_.boundary().findPatchID(patchName_);
    if (patchI == -1)
    {
        FatalError<< "patchName " << patchName_ << " not found in dict"
            << abort(FatalError);
    }

    // const volVectorField& field = mesh_.lookupObject<volVectorField>(fieldName_);
    // fvPatchField<vector> patchField = field.boundaryField()[patchI];

    // Set user defined or default axial profile
    scalarField axialLoc(0);
    if (dict.found("initialAxialLocations"))
    {
        axialLoc = dict.get<scalarField>("initialAxialLocations");
        profileData_ = dict.get<scalarField>("initialData");
    }
    else
    {
        axialLoc.append(0.0);
        profileData_.append
        (
            dict.get<scalar>("initialData")
        );
    }

    // Stringify the intitial profile data from the user
    word initAxialLoc("");
    word initProfileData("");
    forAll(axialLoc, locI)
    {
        initAxialLoc += std::to_string(axialLoc[locI]) + " ";
        initProfileData += std::to_string(profileData_[locI]) + " ";
    }

    // Communicating with the FMU
    commDataLayer& data = commDataLayer::New(time_);

    // Store initial value
    data.storeObj
    (
        initAxialLoc,
        axialLocationsNameFromFMU_,
        commDataLayer::causality::in
    );
    data.storeObj
    (
        initProfileData,
        outputName_,
        commDataLayer::causality::out
    );

    return false;
}

bool Foam::externalIOObject::extAxialPatch::execute()
{
    if (fieldType_ == "vector")
    {
        return executeVectorField();
    }
    else if (fieldType_ == "scalar")
    {
        return executeScalarField();
    }
    else
    {
        FatalErrorInFunction<< "Available fieldType: (vector, scalar)"
            << exit(FatalError);
    }
    return(false);
}


bool Foam::externalIOObject::extAxialPatch::executeScalarField()
{
    // Get commDataLayer
    commDataLayer& data = commDataLayer::New(time_);

    bool isNewStep = false;
    if (relax_ != 1.0 && time_.timeIndex() != time_.startTimeIndex())
    {
        label& isNewStepTemp = data.getObj<label>
        (
            "new_step", commDataLayer::causality::in
        );
        isNewStep = isNewStepTemp == 1 ? true : false;

        Info<< "Relaxing extAxialPatch" << endl;
    }

    // Extract value reference to FMI
    word& out = data.getObj<word>
    (
        outputName_,
        commDataLayer::causality::out
    );

    if (!mesh_.foundObject<volScalarField>(fieldName_))
    {
        Info<< fieldName_ << " not found" << endl;
        return(false);
    }

    // Extract patch field
    const volScalarField& field = mesh_.lookupObject<volScalarField>(fieldName_);
    label patchI = mesh_.boundary().findPatchID(patchName_);
    fvPatchField<scalar> patchField = field.boundaryField()[patchI];


    // Get axial locations from the FMU
    const word axialLocFromFMU = data.getObj<word>
    (
        axialLocationsNameFromFMU_,
        commDataLayer::causality::in
    );

    // Destringify axial locations from FMU
    word value;
    scalarField axialLoc(0);
    std::stringstream ss(axialLocFromFMU);
    while (getline(ss, value, ' '))
    {
        axialLoc.append(std::stod(value));
    }

    const int axialLocSize(axialLoc.size());

    // Reset axial profile in case mismatch between axial loc and profile data
    if (axialLocSize != profileData_.size())
    {
        const scalar initialValue(profileData_[0]);
        WarningIn("Foam::externalIOObject::extAxialPatch::execute()")
            << nl
            << "    Reset axial profile to the initial value " << initialValue
            << nl
            << endl;
        profileData_.clear();
        forAll(axialLoc, locI)
        {
            profileData_.append(initialValue);
        }
    }

    // Compute distance between sampling points
    // scalarList axialLocDistances(0);
    // for (int locI = 0; locI < axialLocSize-1; locI++)
    // {
    //     axialLocDistances.append(axialLoc[locI+1] - axialLoc[locI]);
    // }

    // Update profile data
    word axialList("");
    scalar zPatch(0), dz(0), z(0), dzMin(GREAT), newValue(0);
    forAll(axialLoc, locI)
    {
        // Closest z position to the patch, to be improved with linear interpolation
        z = axialLoc[locI];
        dzMin = GREAT;
        newValue = 0.0;
        forAll(patchField.patch().Cf(), patchI)
        {
            zPatch = patchField.patch().Cf()[patchI] & pinDirection_;
            dz = mag(z - zPatch);

            if (dz < dzMin)
            {
                dzMin = dz;
                newValue = scale_ * patchField[patchI];
            }
        }

        // Relax the result
        if (isNewStep) {
            profileData_[locI] = newValue;
        } else {
            profileData_[locI] = relax_*newValue + (1.0-relax_)*profileData_[locI];
        }

        // Append in the stringified vector list
        axialList += std::to_string(profileData_[locI]) + " ";
    }

    // Set value in the FMI port
    out = axialList;

    //- Write FMI state
    if (mesh_.time().writeTime())
    {
        fmiState_.set("axialLocations", axialLoc);
        fmiState_.set("values", profileData_);
    }

    return false;
}

bool Foam::externalIOObject::extAxialPatch::executeVectorField()
{
    // Get commDataLayer
    commDataLayer& data = commDataLayer::New(time_);

    bool isNewStep = false;
    if (relax_ != 1.0 && time_.timeIndex() != time_.startTimeIndex())
    {
        label& isNewStepTemp = data.getObj<label>
        (
            "new_step", commDataLayer::causality::in
        );
        isNewStep = isNewStepTemp == 1 ? true : false;

        Info<< "Relaxing extAxialPatch" << endl;
    }

    // Extract value reference to FMI
    word& out = data.getObj<word>
    (
        outputName_,
        commDataLayer::causality::out
    );

    // Extract patch field
    const volVectorField& field = mesh_.lookupObject<volVectorField>(fieldName_);
    label patchI = mesh_.boundary().findPatchID(patchName_);
    fvPatchField<vector> patchField = field.boundaryField()[patchI];
    vectorField NrmlVctr = patchField.patch().nf();

    scalarField dotProduct = patchField & NrmlVctr;


    // Get axial locations from the FMU
    const word axialLocFromFMU = data.getObj<word>
    (
        axialLocationsNameFromFMU_,
        commDataLayer::causality::in
    );

    // Destringify axial locations from FMU
    word value;
    scalarField axialLoc(0);
    std::stringstream ss(axialLocFromFMU);
    while (getline(ss, value, ' '))
    {
        axialLoc.append(std::stod(value));
    }

    const int axialLocSize(axialLoc.size());

    // Reset axial profile in case mismatch between axial loc and profile data
    if (axialLocSize != profileData_.size())
    {
        const scalar initialValue(profileData_[0]);
        WarningIn("Foam::externalIOObject::extAxialPatch::execute()")
            << nl
            << "    Reset axial profile to the initial value " << initialValue
            << nl
            << endl;
        profileData_.clear();
        forAll(axialLoc, locI)
        {
            profileData_.append(initialValue);
        }
    }

    // Compute distance between sampling points
    scalarList axialLocDistances(0);
    for (int locI = 0; locI < axialLocSize-1; locI++)
    {
        axialLocDistances.append(axialLoc[locI+1] - axialLoc[locI]);
    }

    // Update profile data
    word axialList("");
    forAll(axialLoc, locI)
    {
        // Closest z position to the patch, to be improved with linear interpolation
        const scalar z(axialLoc[locI]);
        scalar dzMin(1e23);
        scalar newValue(0);
        forAll(patchField.patch().Cf(), i)
        {
            const scalar zPatch(patchField.patch().Cf()[i] & pinDirection_);

            const scalar dz(mag(z - zPatch));

            if (dz < dzMin)
            {
                dzMin = dz;
                newValue = scale_ * mag(dotProduct[i]);
            }
        }

        // Relax the result
        if (isNewStep) {
            profileData_[locI] = newValue;
        } else {
            profileData_[locI] = relax_*newValue + (1.0-relax_)*profileData_[locI];
        }

        // Append in the stringified vector list
        axialList += std::to_string(profileData_[locI]) + " ";
    }

    // Set value in the FMI port
    out = axialList;

    return false;
}

bool Foam::externalIOObject::extAxialPatch::write()
{
    return false;
}


#endif

// ************************************************************************* //
