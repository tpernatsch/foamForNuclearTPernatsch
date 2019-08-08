/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/
#include "subscaleFuel.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fvmSup.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(subscaleFuel, 0);
    defineRunTimeSelectionTable(subscaleFuel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subscaleFuel::subscaleFuel
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "fuelProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    fuelZoneNumber_(PtrList<entry>(IOdictionary::lookup("zones")).size()),
    fuelZoneID_(fuelZoneNumber_),
    Tf_(fuelZoneNumber_),
    Tc_(fuelZoneNumber_),
    fuelSubMeshSize_(fuelZoneNumber_),
    cladSubMeshSize_(fuelZoneNumber_)

{


    PtrList<entry> entries(IOdictionary::lookup("zones"));

    //create 2 fieldFields to initialize Tf_ and Tc_ if not present in the current time in the case folder

    //- Initial solution fields - fuel
    PtrList<FieldField<Field, scalar> > Tf0(fuelZoneNumber_);

    //- Initial  solution fields - clad
    PtrList<FieldField<Field, scalar> > Tc0(fuelZoneNumber_);

    forAll(fuelZoneID_,zoneI)
    {

        const word& name = entries[zoneI].keyword();
        dictionary& dict = entries[zoneI].dict();
        label zoneId = mesh_.cellZones().findZoneID(name);

        fuelZoneID_.set(zoneI, new label(zoneId));

        fuelSubMeshSize_[zoneI] = dict.lookupOrDefault("fuelSubMeshSize",6);
        cladSubMeshSize_[zoneI] = dict.lookupOrDefault("cladSubMeshSize",3);

        Tf0.set(zoneI,new FieldField<Field, scalar>(mesh_.cellZones()[zoneId].size()));
        Tc0.set(zoneI,new FieldField<Field, scalar>(mesh_.cellZones()[zoneId].size()));

        forAll(Tf0[zoneI],cellI)
        {
            Tf0[zoneI].set(cellI, new Field<scalar>(fuelSubMeshSize_[zoneI], dict.lookupOrDefault("Tf0",300.0)));
            Tc0[zoneI].set(cellI, new Field<scalar>(cladSubMeshSize_[zoneI], dict.lookupOrDefault("Tc0",300.0)));
        }

        Tf_.set
        (
            zoneI,
            new IOFieldField<Field,scalar>
            (
                IOobject
                (
                    "Tf_"+mesh_.cellZones()[zoneId].name(),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                Tf0[zoneI]
            )

        );
        Tf_[zoneI].oldTime(); //initialize  Tf0_

        Tc_.set
        (
            zoneI,
            new IOFieldField<Field,scalar>
            (
                IOobject
                (
                    "Tc_"+mesh_.cellZones()[zoneId].name(),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                Tc0[zoneI]
            )

        );
        Tc_[zoneI].oldTime(); //initialize  Tc0_

        if(!Tc_[zoneI].empty() && !Tf_[zoneI].empty())//when running in parallel, you may eg have a list of empty elements if the processor is operating in a region without fuel
        {

            if ((Tc_[zoneI][0].size() !=   cladSubMeshSize_[zoneI]) || (Tf_[zoneI][0].size() !=   fuelSubMeshSize_[zoneI]))
            {
                FatalIOErrorIn("subscaleFuel::"
                    "subscaleFuel(...)", *this)
                    << "subscale mesh size specified in the fuelProperties dict not the same as the one read from the case folder "
                    << abort(FatalIOError);
            }

        }
    }

}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::subscaleFuel> Foam::subscaleFuel::New
(
    const fvMesh& mesh
)
{
    word modelName;

    {
        IOdictionary dict
        (
            IOobject
            (
                "fuelProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dict.lookup("model") >> modelName;
    }

    Info<< "Selecting fuel model type " << modelName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "subscaleFuel::New(const fvMesh&) "
        )   << "Unknown fuel model " << modelName
            << endl << endl
            << "Valid models types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<subscaleFuel>
    (
        cstrIter()(mesh) //here is a crucial point: if a fuel model exists, the functions return its contructor
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::subscaleFuel::~subscaleFuel()
{}


