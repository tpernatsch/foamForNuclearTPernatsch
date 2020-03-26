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

#include "pointKineticNeutronics.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointKineticNeutronics, 0);

    addToRunTimeSelectionTable
    (
        neutronics,
        pointKineticNeutronics,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointKineticNeutronics::pointKineticNeutronics
(
    fvMesh& mesh
)
:
    neutronics(mesh),
    nuclearData_
    (
        IOobject
        (
            "nuclearData",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    power_
    (
        reactorState_.get<scalar>
        (
            "pTarget"
        )
    ),
    externalReactivity_
    (
        reactorState_.lookupOrDefault<scalar>
        (
            "externalReactivity",
            0.0
        )
    ),
    totalReactivity_
    (
        externalReactivity_
    ),
    promptGenerationTime_
    (
        nuclearData_.get<scalar>
        (
            "promptGenerationTime"
        )
    ),
    beta_(0.0),
    betas_
    (
        nuclearData_.get<scalarList>("Beta") // ...
    ),
    lambdas_
    (
        nuclearData_.get<scalarList>("lambda") // ... ...
    ),
    precursorPowers_
    (
        reactorState_.lookupOrDefault<scalarList>
        (
            "precursorPowers",
            scalarList(betas_.size(), 0.0)
        )
    ),
    delayedGroups_(betas_.size()),
    timeIndex_(mesh.time().timeIndex()),
    powerOld_(power_),
    precursorPowersOld_(precursorPowers_),
    coeffTFuel_
    (
        nuclearData_.get<scalar>("feedbackCoeffTFuel")
    ),
    coeffTClad_
    (
        nuclearData_.get<scalar>("feedbackCoeffTClad")
    ),
    coeffTCool_
    (
        nuclearData_.get<scalar>("feedbackCoeffTCool")
    ),
    coeffRhoCool_
    (
        nuclearData_.get<scalar>("feedbackCoeffRhoCool")
    ),
    coeffTStruct_
    (
        nuclearData_.get<scalar>("feedbackCoeffTStruct")
    ),
    TFuel_
    (
        IOobject
        (
            "pointKinetics.TFuel",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TClad_
    (
        IOobject
        (
            "TClad",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TCool_
    (
        IOobject
        (
            "TCool",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    rhoCool_
    (
        IOobject
        (
            "rhoCool",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimDensity, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TStruct_
    (
        IOobject
        (
            "TStruct",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    oneGroupFlux_
    (
        IOobject
        (
            "oneGroupFlux",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless/dimArea/dimTime, 1),
        zeroGradientFvPatchScalarField::typeName
    ),
    energyGroups_(0),
    fluxes_(0),
    precursors_(0),
    liquidFuelPrecursorPowers_(0),
    fastNeutrons_
    (
        nuclearData_.lookupOrDefault<bool>("fastNeutrons", false)
    ),
    coeffFastDoppler_
    (
        (fastNeutrons_) ?
        nuclearData_.get<scalar>("feedbackCoeffFastDoppler") :
        0.0
    ),
    fuelFeedbackCellField_
    (
        IOobject
        (
            "pointKinetics.fuelFeedbackCellField",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    ),
    coolFeedbackCellField_
    (
        IOobject
        (
            "pointKinetics.coolFeedbackCellField",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    ),
    structFeedbackCellField_
    (
        IOobject
        (
            "pointKinetics.structFeedbackCellField",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    )
{
    //- Cannot work in eigenvalue mode for obvious reasons, it makes no sense
    if (eigenvalueNeutronics_)
    {
        FatalErrorInFunction
            << "pointKinetics model incompatible with eigenvalueNeutronics"
            << exit(FatalError);
    }

    if (liquidFuel_)
    {
        FatalErrorInFunction
            << "pointKinetics model does not currently support liquidFuel"
            << exit(FatalError);
    }

    //- Check that provided power is > 0
    if (power_ <= 0)
    {
        FatalErrorInFunction
            << "Set pTarget > 0 !"
            << exit(FatalError);
    }

    //- Check that provided promptGenerationTime is > 0
    if (promptGenerationTime_ <= 0)
    {
        FatalErrorInFunction
            << "Set promptGenerationTime > 0 !"
            << exit(FatalError);
    }

    //- Create fluxes_ PtrList. If pre-existing flux fields exist in the start
    //  time folder, read those and determine the number of energy groups. If
    //  they do not exist, read the number of energy groups from nuclearData
    //  and create the corresponding number of flux fields, scaled so to add
    //  up to oneGroupFlux (if it does not exist, defaults to 1/m2/s). If
    //  the energyGroups keyword is not found, deafults to only one energy
    //  group. Again, most of this has nothing to do with the pointKinetics
    //  model iteself (there is no energy dependence), but it is used to
    //  scale all the fluxes accordingly to pointKinetic results
    while (true)
    {
        word fluxName = "flux"+Foam::name(energyGroups_);

        IOobject fluxHeader
        (
            fluxName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        if (fluxHeader.typeHeaderOk<volScalarField>(true))
        {
            fluxes_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        fluxName,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
            energyGroups_++;
        }
        else
        {
            break;
        }
    }
    if (fluxes_.size() == 0)
    {
        energyGroups_ = 
            nuclearData_.lookupOrDefault<scalar>("energyGroups", 1);

        for (int i = 0; i < energyGroups_; i++)
        {
            fluxes_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        "flux"+Foam::name(i),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        (
                            energyGroups_ == 1 ?
                            IOobject::NO_WRITE :
                            IOobject::AUTO_WRITE
                        )
                    ),
                    (1.0/energyGroups_)*oneGroupFlux_
                )
            );
            fluxes_[i].correctBoundaryConditions();
        }
    }
    else
    {
        oneGroupFlux_ *= 0.0;
        forAllConstIter
        (
            PtrList<volScalarField>,
            fluxes_,
            iter
        )
        {
            const volScalarField& fluxi(iter());
            oneGroupFlux_ += fluxi;
        }
        oneGroupFlux_.correctBoundaryConditions();
    }

    //- Read real precursors if present, they only get re-scaled by 
    //  pointKinetic results
    label i = 0;
    while (true)
    {
        word precName = "prec"+Foam::name(i);

        IOobject precHeader
        (
            precName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        if (precHeader.typeHeaderOk<volScalarField>(true))
        {
            precursors_.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        precName,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
            i++;
        }
        else
        {
            break;
        }
    }

    //- Set precursorPowers so that, if they are not found in reactorState, the
    //  system starts from a steady-state
    if (!reactorState_.found("precursorPowers"))
    {
        for (int i = 0; i < delayedGroups_; i++)
        {
            precursorPowers_[i] = 
                (betas_[i]*power_)/(lambdas_[i]*promptGenerationTime_);
        }
    }

    //- Compute total effective delayed neutron fraction
    forAll(betas_, i)
    {
        beta_ += betas_[i];
    }

    //- Fill in feedback marker fields - Fuel
    wordList fuelFeedbackZones
    (
        nuclearData_.lookupOrDefault("fuelFeedbackZones", wordList())
    );
    if (fuelFeedbackZones.size() == 0)
    {
        forAll(fuelFeedbackCellField_, i)
        {
            fuelFeedbackCellField_[i] = 1.0;
        }
    }
    else
    {
        forAll(fuelFeedbackZones, i)
        {
            word zoneName(fuelFeedbackZones[i]);
            labelList zoneCells(mesh_.cellZones()[zoneName]);
            forAll(zoneCells, j)
            {
                label celli(zoneCells[j]);
                fuelFeedbackCellField_[celli] = 1.0;
            }
        }
    }

    //- Coolant
    wordList coolFeedbackZones
    (
        nuclearData_.lookupOrDefault("coolFeedbackZones", wordList())
    );
    if (coolFeedbackZones.size() == 0)
    {
        forAll(coolFeedbackCellField_, i)
        {
            coolFeedbackCellField_[i] = 1.0;
        }
    }
    else
    {
        forAll(coolFeedbackZones, i)
        {
            word zoneName(coolFeedbackZones[i]);
            labelList zoneCells(mesh_.cellZones()[zoneName]);
            forAll(zoneCells, j)
            {
                label celli(zoneCells[j]);
                coolFeedbackCellField_[celli] = 1.0;
            }
        }
    }

    //- Structures
    wordList structFeedbackZones
    (
        nuclearData_.lookupOrDefault("structFeedbackZones", wordList())
    );
    if (structFeedbackZones.size() == 0)
    {
        forAll(structFeedbackCellField_, i)
        {
            structFeedbackCellField_[i] = 1.0;
        }
    }
    else
    {
        forAll(structFeedbackZones, i)
        {
            word zoneName(structFeedbackZones[i]);
            labelList zoneCells(mesh_.cellZones()[zoneName]);
            forAll(zoneCells, j)
            {
                label celli(zoneCells[j]);
                structFeedbackCellField_[celli] = 1.0;
            }
        }
    }

    //- Some notes on modelling choices, for clarity
    Info<< "The pointKinetics neutronics model currently computes average "
        << "perturbed values for feedback fields (T fuel, cladding, etc.) "
        << "by weighting those over oneGroupFlux. This is not technically "
        << "correct as they should be weighted also via the adjoint flux."
        << "This will be addressed in future updates, but given that the "
        << "adjoint flux is equal to the flux if dealing with only one energy "
        << "group, this is deemed fine for now." << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::pointKineticNeutronics::correct
(
    scalar& residual,
    label couplingIter
) 
{
    #include "solvePointKinetics.H"
}

void Foam::pointKineticNeutronics::getCouplingFieldRefs
(
    const objectRegistry& src,
    const meshToMesh& neutroToFluid
)
{
    //- Field names must reflect those defined in createCouplingFields.H
    TFuelOrig_ = 
        src.findObject<volScalarField>("bafflelessTFuelAv");
    TCladOrig_ = 
        src.findObject<volScalarField>("bafflelessTCladAv");
    TCoolOrig_ = 
        src.findObject<volScalarField>("bafflelessTCool");
    rhoCoolOrig_ = 
        src.findObject<volScalarField>("bafflelessRhoCool");
    
    //- Project thermalHydraulic volFuelPower onto the neutronic one to
    //  initialize it if the latter does not exist
    IOobject volFuelPowerHeader
    (
        "volFuelPower",
        mesh_.time().timeName(),
        mesh_.time(),
        IOobject::NO_READ
    );

    if (!volFuelPowerHeader.typeHeaderOk<volScalarField>(true))
    {
        volFuelPowerOrig_ = 
            src.findObject<volScalarField>("bafflelessVolFuelPower");
        neutroToFluid.mapTgtToSrc
            (
                *volFuelPowerOrig_, 
                plusEqOp<scalar>(), 
                volFuelPower_
            );
        volFuelPower_.correctBoundaryConditions(); 
    }

    //- These have no use now, will be needed by the liquid
    /*
    if (liquidFuel_)
    {
        UOrig_ = 
            src.findObject<volVectorField>("bafflelessU");
        alphaOrig_ = 
            src.findObject<volScalarField>("bafflelessAlpha");
        alphatOrig_ = 
            src.findObject<volScalarField>("bafflelessAlphat");
        muOrig_ =
            src.findObject<volScalarField>("bafflelessMu");
    }
    else
    {
        UOrig_ = nullptr;
        alphaOrig_ = nullptr;
        alphatOrig_ = nullptr;
        muOrig_ = nullptr;
    }
    */

    //- The rest of this function is for initializing the reference values of
    //  the feedback parameters. If they are found in the dictionary, use
    //  those, otherwise compute them from the coupling fields (thus
    //  assuming that the simulation starts from a steady state)
    
    this->interpolateCouplingFields(neutroToFluid);

    #include "computeFeedbackFieldValues.H"

    TFuelRef_ = 
        max
        (
            nuclearData_.lookupOrDefault<scalar>("TFuelRef", TFuelValue),
            SMALL
        );
    
    TCladRef_ =
        nuclearData_.lookupOrDefault<scalar>("TCladdingRef", TCladValue);

    TCoolRef_ =
        nuclearData_.lookupOrDefault<scalar>("TCoolantRef", TCoolValue);

    rhoCoolRef_ =
        nuclearData_.lookupOrDefault<scalar>("rhoCoolantRef", rhoCoolValue);

    TStructRef_ =
        nuclearData_.lookupOrDefault<scalar>("TStructuresRef", TStructValue);

    #include "correctReactivity.H"

    Info << endl << "pointKinetics (initial conditions): " << endl;
    Info << "    power = " << power_ << " W" << endl;
    Info << "    totalReactivity   = " << (totalReactivity_*1e5)<< " pcm" 
                                                                << endl;
    Info << "    -> Doppler (fast) = " << (DopplerReactivity*1e5)<< " pcm" 
                                                                << endl;
    Info << "    -> TFuel          = " << (TFuelReactivity*1e5) << " pcm" 
                                                                << endl;
    Info << "    -> TClad          = " << (TCladReactivity*1e5) << " pcm" 
                                                                << endl;
    Info << "    -> TCool          = " << (TCoolReactivity*1e5) << " pcm" 
                                                                << endl;
    Info << "    -> rhoCool        = " << (rhoCoolReactivity*1e5) << " pcm" 
                                                                << endl;
    Info << "    -> TStruct        = " << (TStructReactivity*1e5) <<" pcm" 
                                                                << endl;
    Info << "    TFuel = " << TFuelValue << " K" << endl;
    Info << "    TClad = " << TCladValue << " K" << endl;
    Info << "    TCool = " << TCoolValue << " K" << endl;
    Info << "    rhoCool = " << rhoCoolValue << " kg/m3" << endl;
    Info << "    TStruct = " << TStructValue << " K" << endl;
}

void Foam::pointKineticNeutronics::interpolateCouplingFields
(
    const meshToMesh& neutroToFluid
)
{
    neutroToFluid.mapTgtToSrc(*TFuelOrig_, plusEqOp<scalar>(), TFuel_);
    neutroToFluid.mapTgtToSrc(*TCladOrig_, plusEqOp<scalar>(), TClad_);
    neutroToFluid.mapTgtToSrc(*TCoolOrig_, plusEqOp<scalar>(), TCool_);
    neutroToFluid.mapTgtToSrc(*rhoCoolOrig_, plusEqOp<scalar>(), rhoCool_);
    
    /* This is for when liquidFuel support will be added
    if (liquidFuel_)
    {
        neutroToFluid.mapTgtToSrc(*UOrig_, plusEqOp<vector>(), UPtr_());
        neutroToFluid.mapTgtToSrc
        (
            *alphaOrig_, 
            plusEqOp<scalar>(), 
            alphaPtr_()
        );
        neutroToFluid.mapTgtToSrc
        (
            *alphatOrig_, 
            plusEqOp<scalar>(), 
            alphatPtr_()
        );
        neutroToFluid.mapTgtToSrc(*muOrig_, plusEqOp<scalar>(), muPtr_());
        phiPtr_() = fvc::flux(UPtr_());
        volScalarField diffCoeffOrig
        (
            (
                *alphatOrig_ 
            +   *muOrig_/xs_.ScNo()
            )/(*rhoCoolOrig_)
        ); 
        neutroToFluid.mapTgtToSrc
        (
            diffCoeffOrig, 
            plusEqOp<scalar>(), 
            diffCoeffPrecPtr_()
        );

        UPtr_().correctBoundaryConditions();
        alphaPtr_().correctBoundaryConditions();
        alphatPtr_().correctBoundaryConditions();
        diffCoeffPrecPtr_().correctBoundaryConditions();
    }
    */

    TFuel_.correctBoundaryConditions();
    TClad_.correctBoundaryConditions();
    TCool_.correctBoundaryConditions();
    rhoCool_.correctBoundaryConditions();
}

// ************************************************************************* //
