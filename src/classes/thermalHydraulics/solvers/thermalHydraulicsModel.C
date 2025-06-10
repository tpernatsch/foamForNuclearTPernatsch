/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2412                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2024 OpenCFD Ltd.          |
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

\*---------------------------------------------------------------------------*/

#include "thermalHydraulicsModel.H"
#include "regimeMapModel.H"
#include "mergeOrSplitBaffles.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(thermalHydraulicsModel, 0);
    // defineRunTimeSelectionTable(thermalHydraulicsModel, thermalHydraulicsModels);
}
}

const Foam::Enum
<
    Foam::solvers::thermalHydraulicsModel::momentumMode
>
Foam::solvers::thermalHydraulicsModel::momentumModeNames_
(
    {
        {
            momentumMode::cellCentered,
            "cellCentered"
        },
        {
            momentumMode::cellCenteredFaceReconstruction,
            "cellCenteredFaceReconstruction"
        },
        {
            momentumMode::faceCentered,
            "faceCentered"
        }
    }
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::thermalHydraulicsModel::thermalHydraulicsModel
(
    const Time& time,
    dynamicFvMesh& mesh,
    fv::options& fvOptions
)
:
    solver(mesh),
    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    runTime_(time),
    mesh_(mesh),
    pimple_(mesh),
    fvOptions_(fvOptions),
    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimVol/dimTime, 0.0)
    ),
    g_
    (
        IOobject
        (
            "g",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    hRef_
    (
        IOobject
        (
            "hRef",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar("hRef", dimLength, 0)
    ),
    gh_
    (
        "gh",
        (g_ & mesh.C()) + mag(g_)*hRef_
    ),
    ghf_
    (
        "ghf",
        (g_ & mesh.Cf()) + mag(g_)*hRef_
    ),
    p_rgh_
    (
        IOobject
        (
            "p_rgh",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        p_rgh_
    ),
    stressTensor_
    (
        IOobject
        (
            "stressTensor",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("", dimPressure, symmTensor::zero)
    ),
    kappaEff_
    (
        IOobject
        (
            "kappaEff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimPower/dimTemperature/dimLength, 0)
    ),
    pMin_
    (
        "pMin",
        p_.dimensions(),
        *this
    ),
    pRefCell_
    (
        this->lookupOrDefault("pRefCell", 0)
    ),
    pRefValue_
    (
        this->lookupOrDefault("pRefValue", 0)
    ),
    forcePRef_
    (
        this->lookupOrDefault<bool>("forcePRef", false)
    ),
    initialFluidMass_("initialFluidMass", dimMass, 0),
    momentumMode_
    (
        momentumModeNames_.get
        (
            pimple_.dict().lookupOrDefault<word>
            (
                "momentumMode",
                "cellCentered"
            )
        )
    ),
    porousInterfaceSharpness_
    (
        pimple_.dict().lookupOrDefault<scalar>
        (
            "porousInterfaceSharpness",
            0.0
        )
    ),
    powerDensityNeutronics_
    (
        IOobject
        (
            "powerDensityNeutronics",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower/dimVol, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    powerDensityNeutronicsToLiquid_
    (
        IOobject
        (
            "powerDensityNeutronicsToLiquid",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower/dimVol, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    originalPoints_(mesh_.points())
{
    -
    setRefCell
    (
        p_,
        p_rgh_,
        *this,
        pRefCell_,
        pRefValue_,
        forcePRef_
    );
    mesh.setFluxRequired(p_rgh_.name());

    //- Construct regime maps, if any
    if (this->found("regimeMapModels"))
    {
        const dictionary& regimeMapModelsDict
        (
            this->subDict("regimeMapModels")
        );
        wordList regimeMapNames = regimeMapModelsDict.toc();
        forAll(regimeMapNames, i)
        {
            const dictionary& regimeMapModelDict
            (
                regimeMapModelsDict.subDict(regimeMapNames[i])
            );
            regimeMapModels_.set
            (
                regimeMapModelDict.dictName(),
                autoPtr<regimeMapModel>
                (
                    regimeMapModel::New
                    (
                        mesh_,
                        regimeMapModelDict
                    )
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::thermalHydraulicsModel::~thermalHydraulicsModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::solvers::thermalHydraulicsModel::correctRegimeMaps()
{
    forAllIter
    (
        regimeMapTable,
        regimeMapModels_,
        iter
    )
    {
        iter()->correct();
    }
}

void Foam::solvers::thermalHydraulicsModel::stop()
{
    FatalErrorInFunction
        << "Terminating execution"
        << exit(FatalError);
}

void Foam::solvers::thermalHydraulicsModel::correctBaffleLessFields()
{
    if (mesh_.time().controlDict().found("removeBaffles"))
    {
        const dictionary& removeBafflesDict = mesh_.time().controlDict().subDict("removeBaffles");
        if (removeBafflesDict.getOrDefault<bool>(mesh_.name(), false))
        {

            const IOdictionary couplingDict
            (
                IOobject
                (
                    "regionsDict",
                    runTime.time().constant(),
                    runTime.db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            );

            //Lookup for the fields that need to be mapped FROM this mesh
            const dictionary mappingDict(couplingDict.subDict("mappings"));
            //Loop on every region that is not this one and look for the fields in "sourceFields"
            const wordList regions(mappingDict.toc());

            forAll(regions, regioni)
            {
                if (regions[regioni] != mesh_.name()) //look for other regions
                {
                    const dictionary regionFromDict(mappingDict.subDict(regions[regioni]));
                    const wordList regionsFrom(regionFromDict.toc());
                    forAll(regionsFrom, regionFromi)
                    {
                        if (regionsFrom[regionFromi] == mesh_.name())
                        {
                            const wordList fieldsList(regionFromDict.subDict(regionsFrom[regionFromi]).get<wordList>("sourceFields")); // list of fields to create
                            fvMesh& baffleLessMesh = const_cast<fvMesh&>(mesh_.time().lookupObject<fvMesh>(mesh_.name()+".baffleLess"));
                            forAll(fieldsList, fieldi)
                            {
                                correctBaffleLessField<scalar>(fieldsList[fieldi], baffleLessMesh);
                                correctBaffleLessField<vector>(fieldsList[fieldi], baffleLessMesh);
                                correctBaffleLessField<tensor>(fieldsList[fieldi], baffleLessMesh);
                                correctBaffleLessField<symmTensor>(fieldsList[fieldi], baffleLessMesh);
                                correctBaffleLessField<sphericalTensor>(fieldsList[fieldi], baffleLessMesh);

                            }
                        }
                    }
                }
            }
        }
    }
}

template<class Type>
void Foam::solvers::thermalHydraulicsModel::correctBaffleLessField(word fieldName, fvMesh& baffleLessMesh)
{

    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    if(mesh_.foundObject<VolFieldType>(fieldName))
    {
        VolFieldType& field = baffleLessMesh.lookupObjectRef<VolFieldType>(fieldName+".baffleLess");
        field.primitiveFieldRef() = mesh_.lookupObject<VolFieldType>(fieldName).primitiveField();
        field.correctBoundaryConditions();
    }
}

void Foam::solvers::thermalHydraulicsModel::deformMesh()
{
    //-Look for the multiRegionDict

    const IOdictionary couplingDict
    (
        IOobject
        (
            "regionsDict",
            runTime.time().constant(),
            runTime.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    if (couplingDict.found("meshDeformation"))
    {

        const dictionary deformDict(couplingDict.subDict("meshDeformation"));

        const wordList regions(deformDict.toc());

        forAll(regions, regioni)
        {
            if (regions[regioni] == mesh_.name())
            {
                const volPointInterpolation& meshPointInterpolation = volPointInterpolation::New(mesh_);

                tmp<pointVectorField> meshPointsDisplacement = meshPointInterpolation.interpolate
                (
                    mesh_.lookupObject<volVectorField>
                    (
                        deformDict.subDict(mesh_.name()).get<word>("displacementField")
                    )
                );

                tmp<pointField> displacedPoints = originalPoints_ + meshPointsDisplacement->internalField();

                mesh_.movePoints(displacedPoints);
            }
        }
    }
}

// void Foam::thermalHydraulicsModel::adjustTimeStep()
// {
//     this->correctCourant();

//     bool adjustTimeStep =
//         runTime_.controlDict().lookupOrDefault("adjustTimeStep", false);

//     scalar maxCo =
//         runTime_.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

//     scalar maxDeltaT =
//         runTime_.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

//     if (adjustTimeStep)
//     {
//         scalar maxDeltaTFact = maxCo/(CoNum_ + SMALL);
//         scalar deltaTFact =
//             min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

//         runTime_.setDeltaT
//         (
//             min
//             (
//                 deltaTFact*runTime_.deltaTValue(),
//                 maxDeltaT
//             )
//         );
//     }
// }





// ************************************************************************* //
