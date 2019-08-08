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

#include "porousMedium.H"
#include "zeroGradientFvPatchFields.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(porousMedium, 0);
    defineRunTimeSelectionTable(porousMedium, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousMedium::porousMedium
(
    const volScalarField& rho,
    const volVectorField& U,
    const rhoThermo& thermo
)
:
    IOdictionary(// NB: porousMedium deruved from IOdictionary: need to construct IOdictionary
        IOobject(// NB: IOdictionary need IOobject as constructor
            "porousMediumProperties",
            U.mesh().time().constant(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(U.mesh()),
    thermo_(thermo),
    rho_(rho),
    U_(U),
    mu_(thermo.mu()),
    gamma_(
        IOobject(
            "gamma",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(word(), dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    volumetricAreaSS_(
        IOobject(
            "volumetricAreaSS",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(word(), dimensionSet(0,-1,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    volumetricAreaFuel_(
        IOobject(
            "volumetricAreaFuel",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(word(), dimensionSet(0,-1,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    volumetricAreaHX_(
        IOobject(
            "volumetricAreaHX",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(word(), dimensionSet(0,-1,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    rotate_(
        IOobject(
            "porousMedium::rotate",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(word(), dimless, tensor::I),
        zeroGradientFvPatchScalarField::typeName
    ),
    externalHeatTransferCoefficient_(
        IOobject(
            "externalHeatTransferCoefficient",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(word(), dimensionSet(1,0,-3,-1,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    externalT_(
        IOobject(
            "externalT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(word(), dimensionSet(0,0,0,1,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    externalVolHeatSource_(
        IOobject(
            "externalVolHeatSource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(word(), dimensionSet(1,-1,-3,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    externalRhoCp_(
        IOobject(
            "externalRhoCp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(word(), dimensionSet(1,-1,-2,-1,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    )
{ }

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::porousMedium> Foam::porousMedium::New
(
    const volScalarField& rho,
    const volVectorField& U,
    const rhoThermo& thermo
)
{
    word modelName;

    {
        IOdictionary dict
        (
            IOobject
            (
                "porousMediumProperties",
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dict.lookup("model") >> modelName;
    }

    Info<< "Selecting porous medium model type " << modelName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "porousMedium::New(const volScalarField&, "
            "const volVectorField&, basicThermo&)"
        )   << "Unknown porous medium model " << modelName
            << endl << endl
            << "Valid models types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<porousMedium>
    (
        cstrIter()(rho, U, thermo) //if a porousMedium model exists, the functions return its contructor
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porousMedium::~porousMedium()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::tmp< Foam::volVectorField >
Foam::porousMedium::explicitMomentumSource() const //totally explicit
{
    //- Friction coefficient
    return (-gamma_*dragCoeff() & U_) + pump();
}


Foam::tmp<Foam::fvVectorMatrix>  
Foam::porousMedium::semiImplicitMomentumSource(volVectorField& U) const  //linear dependence on U left implicit
{

    //- Friction coefficient
   return -gamma_*(fvm::Sp(1.0/3.0*tr(dragCoeff()) / gamma_, U) + (dev(dragCoeff()) & U / gamma_)) + pump();
}

void Foam::porousMedium::correct()
{}

// ************************************************************************* //
