/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "fluid.H"
#include "fvCFD.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "fvcSnGrad.H"
#include "fvcCurl.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedValueFvsPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::fluid::stateOfMatter
>
Foam::fluid::stateOfMatterNames_
(
    {
        { 
            stateOfMatter::undetermined, 
            "undetermined" 
        },
        { 
            stateOfMatter::liquid, 
            "liquid" 
        },
        { 
            stateOfMatter::gas, 
            "gas" 
        }
    }
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluid::fluid
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& phaseName,
    bool readIfPresentAndWrite
)
:
    phaseBase
    ( 
        dict,
        mesh,
        phaseName,
        readIfPresentAndWrite,
        readIfPresentAndWrite
    ),
    stateOfMatter_
    (
        stateOfMatterNames_.get
        (
            dict_.lookupOrDefault<word>
            (
                "stateOfMatter", 
                "undetermined"
            )
        )
    ),
    U_
    (
        IOobject
        (
            IOobject::groupName("U", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("", dimVelocity, vector(0,0,0)),
        slipFvPatchVectorField::typeName
    ),
    magU_
    (
        IOobject
        (
            IOobject::groupName("magU", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimVelocity, scalar(0))
    ),
    alphaPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            (
                (
                    mesh.time().controlDict().lookupOrDefault<bool>
                    (
                        "writeRestartFields", 
                        true
                    )
                ) ?
                IOobject::AUTO_WRITE :
                IOobject::NO_WRITE
            )
        ),
        mesh,
        dimensionedScalar("", dimVol/dimTime, 0)
    ),
    alphaRhoPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPhi", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            (
                (
                    mesh.time().controlDict().lookupOrDefault<bool>
                    (
                        "writeRestartFields", 
                        true
                    )
                ) ?
                IOobject::AUTO_WRITE :
                IOobject::NO_WRITE
            )
        ),
        mesh,
        dimensionedScalar("", dimMass/dimTime, 0)
    ),
    dgdt_
    (
        IOobject
        (
            IOobject::groupName("dgdt", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            (
                (
                    mesh.time().controlDict().lookupOrDefault<bool>
                    (
                        "writeRestartFields", 
                        true
                    )
                ) ?
                IOobject::AUTO_WRITE :
                IOobject::NO_WRITE
            )
        ),
        mesh,
        dimensionedScalar("", dimless/dimTime, 0)
    ),
    staticQuality_
    (
        IOobject
        (
            IOobject::groupName("staticQuality", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            (
                (
                    mesh.time().controlDict().lookupOrDefault<bool>
                    (
                        "writeAllFields", 
                        false
                    )
                ) ?
                IOobject::AUTO_WRITE :
                IOobject::NO_WRITE
            )
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    ),
    flowQuality_
    (
        IOobject
        (
            IOobject::groupName("flowQuality", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            (
                (
                    mesh.time().controlDict().lookupOrDefault<bool>
                    (
                        "writeAllFields", 
                        false
                    )
                ) ?
                IOobject::AUTO_WRITE :
                IOobject::NO_WRITE
            )
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    ),
    normalized_
    (
        IOobject
        (
            IOobject::groupName("normalized.alpha", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    ),
    Pr_
    (
        IOobject
        (
            IOobject::groupName("Pr", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    ),
    contErr_
    (
        IOobject
        (
            IOobject::groupName("contErr", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            (
                (
                    mesh.time().controlDict().lookupOrDefault<bool>
                    (
                        "writeContinuityErrors", 
                        false
                    )
                ) ?
                IOobject::AUTO_WRITE :
                IOobject::NO_WRITE
            )
        ),
        mesh,
        dimensionedScalar("", dimDensity/dimTime, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Dh_
    (
        IOobject
        (
            IOobject::groupName("Dh", this->name()),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimLength, SMALL),
        zeroGradientFvPatchScalarField::typeName
    ),
    dispersion_(mesh.C().size(), scalar(0.0)),
    thermoResidualAlpha_
    (
        dimensionedScalar::lookupOrDefault
        (
            "thermoResidualAlpha",
            dict_,
            dimless, 
            0.0
        )
    )
{
    if (phaseName != "")
    {
        Info << endl << "Constructing fluid: " << this->name() << endl;
    }
    else
    {
        Info << endl << "Constructing fluid" << endl;
    }

    //- Construct thermodynamics package
    thermo_.reset(rhoThermo::New(mesh, this->name()));

    mesh.setFluxRequired(this->name());

    //- Set initial cellZone phase fractions, if present
    if (dict_.found("initialAlphas"))
    {
        const dictionary& alpha0s
        (
            dict_.subDict("initialAlphas")
        );

        forAllConstIter
        (
            dictionary,
            alpha0s,
            alpha0Iter
        )
        {
            labelList cells;
            word zoneName(alpha0Iter->keyword());
            scalar alpha0(alpha0s.get<scalar>(zoneName));
            forAllConstIter
            (
                labelList,
                mesh.cellZones()[zoneName],
                cIter
            )
            {
                cells.append(*cIter);
            }
            forAll(cells, i)
            {
                (*this)[cells[i]] = alpha0;
            }
        }

        this->correctBoundaryConditions();
    }

    thermo_->validate(phaseName, "h", "e");

    //- Init magU and Pr
    magU_ = mag(U_);
    magU_.correctBoundaryConditions();
    Pr_ = thermo_->Cp()*thermo_->mu()/thermo_->kappa();
    Pr_.correctBoundaryConditions();

    //- The rest of the constructor is only for correctly setting the
    //  boundary conditions of phi
    const word phiName = IOobject::groupName("phi", this->name());

    IOobject phiHeader
    (
        phiName,
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ
    );

    if (phiHeader.typeHeaderOk<surfaceScalarField>(true))
    {
        Info<< "Reading face flux field " << phiName << endl;

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }
    else
    {
        Info<< "Calculating face flux field " << phiName << endl;

        wordList phiTypes
        (
            U_.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        forAll(U_.boundaryField(), i)
        {
            if
            (
                isA<fixedValueFvPatchVectorField>(U_.boundaryField()[i])
             || isA<slipFvPatchVectorField>(U_.boundaryField()[i])
             || isA<partialSlipFvPatchVectorField>(U_.boundaryField()[i])
            )
            {
                phiTypes[i] = fixedValueFvsPatchScalarField::typeName;
            }
        }

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    (
                        (
                            mesh.time().controlDict().lookupOrDefault<bool>
                            (
                                "writeRestartFields", 
                                true
                            )
                        ) ?
                        IOobject::AUTO_WRITE :
                        IOobject::NO_WRITE
                    )
                ),
                fvc::flux(U_),
                phiTypes
            )
        );
    }
    
    //- What about alphaPhi, alphaRhoPhi? Well, these depend on the 
    //  phase fraction (unlike phi) but at this step there have been no
    //  phase fractions normalizations (which can only be done by the main
    //  solver). Thus, rather than tentatively set the fields twice (here,
    //  maybe wrong, and then in the main to correct for potential phase
    //  fraction normalization issues), this needs to be handled by the
    //  main solver via the initAlphaPhis function. Needless to say, if
    //  alphaPhi and alphaRhoPhi are found on disk, those are read and that's
    //  the end of it
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluid::~fluid()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluid::initAlphaPhis()
{
    IOobject alphaPhiHeader
    (
        alphaPhi_.name(),
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ
    );
    IOobject alphaRhoPhiHeader
    (
        alphaRhoPhi_.name(),
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ
    );

    if (!alphaPhiHeader.typeHeaderOk<surfaceScalarField>(true))
    {
        alphaPhi_ = fvc::interpolate(*this)*phiPtr_();
    }
    if (!alphaRhoPhiHeader.typeHeaderOk<surfaceScalarField>(true))
    {
        alphaRhoPhi_ = fvc::interpolate(thermo_->rho())*alphaPhi_;
    }
}

void Foam::fluid::constructTurbulenceModel()
{
    turbulence_ =
        phaseCompressibleTurbulenceModel::New
        (
            *this,
            thermo_->rho(),
            U_,
            alphaRhoPhi_,
            phi(),
            thermo_
        ); 
}


// ************************************************************************* //