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

#include "thermoMechanics.H"
#include "zeroGradientFvPatchFields.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermoMechanics, 0);
    defineRunTimeSelectionTable(thermoMechanics, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermoMechanics::thermoMechanics
(
    fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "thermoMechanicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    initialResidual_(1.0),
    meshDisp_
    (
        IOobject
        (
            "meshDisp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("", dimensionSet(0,1,0,0,0,0,0), vector(0,0,0)),
        zeroGradientFvPatchScalarField::typeName
    )
{
Info << "here" << endl;

}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::thermoMechanics> Foam::thermoMechanics::New
(
    fvMesh& mesh
)
{
    word modelName;

    {
        IOdictionary dict
        (
            IOobject
            (
                "thermoMechanicalProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dict.lookup("model") >> modelName;
    }

    Info<< "Selecting thermoMechanics model type " << modelName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "thermoMechanics::New(const volScalarField&, "
            "const volVectorField&, basicThermo&)"
        )   << "Unknown thermoMechanics model " << modelName
            << endl << endl
            << "Valid models types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<thermoMechanics>
    (
        cstrIter()(mesh) //if a thermoMechanics model exists, the functions return its contructor
    );

}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

	bool Foam::thermoMechanics::writeData(Ostream& os) const
	{

	    return os.good();
	}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermoMechanics::~thermoMechanics()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //




// ************************************************************************* //
