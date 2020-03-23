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
    coeffTCladding_
    (
        nuclearData_.get<scalar>("feedbackCoeffTCladding")
    ),
    coeffTCoolant_
    (
        nuclearData_.get<scalar>("feedbackCoeffTCoolant")
    ),
    coeffRhoCoolant_
    (
        nuclearData_.get<scalar>("feedbackCoeffRhoCoolant")
    ),
    coeffTStructures_
    (
        nuclearData_.get<scalar>("feedbackCoeffTStructures")
    ),
    TFuelRef_
    (
        nuclearData_.lookupOrDefault<scalar>("TFuelRef", 0)
    ),
    TCladdingRef_
    (
        nuclearData_.lookupOrDefault<scalar>("TCladdingRef", 0)
    ),
    TCoolantRef_
    (
        nuclearData_.lookupOrDefault<scalar>("TCoolantRef", 0)
    ),
    rhoCoolantRef_
    (
        nuclearData_.lookupOrDefault<scalar>("rhoCoolantRef", 0)
    ),
    TStructuresRef_
    (
        nuclearData_.lookupOrDefault<scalar>("TStructuresRef", 0)
    ),
    resetTFuelRef_
    (
        nuclearData_.found("TFuelRef") ? false : true
    ),
    resetTCladdingRef_
    (
        nuclearData_.found("TCladdingRef") ? false : true
    ),
    resetTCoolantRef_
    (
        nuclearData_.found("TCoolantRef") ? false : true
    ),
    resetRhoCoolantRef_
    (
        nuclearData_.found("rhoCoolantRef") ? false : true
    ),
    resetTStructuresRef_
    (
        nuclearData_.found("TStructuresRef") ? false : true
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
    TCladding_
    (
        IOobject
        (
            "pointKinetics.TCladding",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TCoolant_
    (
        IOobject
        (
            "pointKinetics.TCoolant",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    rhoCoolant_
    (
        IOobject
        (
            "pointKinetics.rhoCoolant",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimDensity, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TStructures_
    (
        IOobject
        (
            "pointKinetics.TStructures",
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
    coolantFeedbackCellField_
    (
        IOobject
        (
            "pointKinetics.coolantFeedbackCellField",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    ),
    structuresFeedbackCellField_
    (
        IOobject
        (
            "pointKinetics.structuresFeedbackCellField",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    ),
    phi_
    (
        IOobject
        (
            "pointKinetics.phi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimVolume/dimTime, 0)
    ),
    flowFaces_(0),
    flowRef_(0),
    flowFraction_(1),
    flowFractionRef_(nuclearData_.get<scalar>("flowFractionRef")),
    sodiumLevelRef_
    (
        nuclearData_.lookupOrDefault<scalar>("sodiumLevelRef", 0.0)
    ),
    GEMReactivityMap_
    (
        nuclearData_.lookupOrDefault<List<Pair<scalar>>>
        (
            "GEMReactivityMap",
            List<Pair<scalar>>()
        )
    ),
    drivelineFeedbackCellField_
    (
        IOobject
        (
            "pointKinetics.drivelineFeedbackCellField",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless, 0)
    ),
    drivelineExpansionCoeff_
    (
        nuclearData_.lookupOrDefault<scalar>("drivelineExpansionCoeff", 0.0)
    ),
    CRReactivityMap_
    (
        nuclearData_.lookupOrDefault<List<Pair<scalar>>>
        (
            "CRReactivityMap",
            List<Pair<scalar>>()
        )
    )
{
    //- Cannot work in eigenvalue mode for obvious reasons, it makes no sense
    if (eigenvalueNeutronics_)
    {
        FatalErrorInFunction
            << "pointKinetics model incompatible with eigenvalueNeutronics!"
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
    wordList coolantFeedbackZones
    (
        nuclearData_.lookupOrDefault("coolantFeedbackZones", wordList())
    );
    if (coolantFeedbackZones.size() == 0)
    {
        forAll(coolantFeedbackCellField_, i)
        {
            coolantFeedbackCellField_[i] = 1.0;
        }
    }
    else
    {
        forAll(coolantFeedbackZones, i)
        {
            word zoneName(coolantFeedbackZones[i]);
            labelList zoneCells(mesh_.cellZones()[zoneName]);
            forAll(zoneCells, j)
            {
                label celli(zoneCells[j]);
                coolantFeedbackCellField_[celli] = 1.0;
            }
        }
    }

    //- Structures
    wordList structuresFeedbackZones
    (
        nuclearData_.lookupOrDefault("structuresFeedbackZones", wordList())
    );
    if (structuresFeedbackZones.size() == 0)
    {
        forAll(structuresFeedbackCellField_, i)
        {
            structuresFeedbackCellField_[i] = 1.0;
        }
    }
    else
    {
        forAll(structuresFeedbackZones, i)
        {
            word zoneName(structuresFeedbackZones[i]);
            labelList zoneCells(mesh_.cellZones()[zoneName]);
            forAll(zoneCells, j)
            {
                label celli(zoneCells[j]);
                structuresFeedbackCellField_[celli] = 1.0;
            }
        }
    }

    //- GEM related
    //- Construct flowFaces_
    wordList flowFaceZoneNames(nuclearData_.get<wordList>("flowFaceZones"));
    forAll(flowFaceZoneNames, i)
    {
        labelList flowFaces(mesh_.faceZones()[flowFaceZoneNames[i]]);
        forAll(flowFaces, j)
        {
            flowFaces_.append(flowFaces[j]);
        }
    }

    //- Check that GEMReactivityMap is indexed by descending sodium level 
    //  values
    for (int i = 0; i < GEMReactivityMap_.size()-1; i++)
    {
        Pair<scalar> X0(GEMReactivityMap_[i]);
        Pair<scalar> X1(GEMReactivityMap_[i+1]);
        if (X0.first() <= X1.first())
        {
            FatalErrorInFunction
            << "GEMReactivityMap should be indexed by descending sodium level"
            << " values!" << exit(FatalError);
        }
    }

    //- Init flowFraction_ to flowFractionRef_
    flowFraction_ = flowFractionRef_;

    //- Driveline related
    wordList drivelineFeedbackZones
    (
        nuclearData_.lookupOrDefault("drivelineFeedbackZones", wordList())
    );
    if (drivelineFeedbackZones.size() == 0)
    {
        forAll(drivelineFeedbackCellField_, i)
        {
            drivelineFeedbackCellField_[i] = 1.0;
        }
    }
    else
    {
        forAll(drivelineFeedbackZones, i)
        {
            word zoneName(drivelineFeedbackZones[i]);
            labelList zoneCells(mesh_.cellZones()[zoneName]);
            forAll(zoneCells, j)
            {
                label celli(zoneCells[j]);
                drivelineFeedbackCellField_[celli] = 1.0;
            }
        }
    }

    //- Add entries to reactivityValues_
    reactivityValues_.insert
    (
        "total",
        0.0
    );
    reactivityValues_.insert
    (
        "external",
        0.0
    );
    reactivityValues_.insert
    (
        "Doppler",
        0.0
    );
    reactivityValues_.insert
    (
        "GEM",
        0.0
    );
    reactivityValues_.insert
    (
        "TFuel",
        0.0
    );
    reactivityValues_.insert
    (
        "TCladding",
        0.0
    );
    reactivityValues_.insert
    (
        "TCoolant",
        0.0
    );
    reactivityValues_.insert
    (
        "TCoolant",
        0.0
    );
    reactivityValues_.insert
    (
        "rhoCoolant",
        0.0
    );
    reactivityValues_.insert
    (
        "TStructures",
        0.0
    );
    reactivityValues_.insert
    (
        "driveline",
        0.0
    );  

    //- Some notes on modelling choices, for clarity
    Info<< "The pointKinetics neutronics model currently computes average "
        << "perturbed values for feedback fields (T fuel, cladding, etc.) "
        << "by weighting those over oneGroupFlux. This is not technically "
        << "correct as they should be weighted also via the adjoint flux."
        << "This will be addressed in future commits, but given that the "
        << "adjoint flux is equal to the flux if dealing with only one energy "
        << "group, this is fine for now." << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::pointKineticNeutronics::correct
(
    const label couplingIter, 
    scalar& residual, 
    const bool& liquidFuel
) 
{
    #include "correctReactivityFeedbacks.H"
    
    if (!liquidFuel)
    {
        #include "solvePointKinetics.H"
    }
    else
    {
        #include "solveLiquidFuelPointKinetics.H"
    }
}

void Foam::pointKineticNeutronics::getFields
(
    const volVectorFieldTable& vvfTable,
    const volScalarFieldTable& vsfTable, 
    const meshToMesh& neutroToFluid
)
{
    //- This if is to compensate for the very debatable implementation
    //  of how volFuelPower is dealt with if the constSubscale fuel power
    //  is provided instead of a volFuelPower (this also required changes
    //  in the constSubscale fuel class)
    if  
    (
            mesh_.time().timeIndex()
        -   mesh_.time().startTimeIndex() == 1 
    )
    {
        neutroToFluid.mapTgtToSrc
        (
            vsfTable["volFuelPower"], 
            plusEqOp<scalar>(), 
            volFuelPower_
        );
        volFuelPower_.correctBoundaryConditions();
    }

    neutroToFluid.mapTgtToSrc(vsfTable["TFuel"], plusEqOp<scalar>(), TFuel_);
    neutroToFluid.mapTgtToSrc(vsfTable["TClad"], plusEqOp<scalar>(), TCladding_);
    neutroToFluid.mapTgtToSrc
    (
        vsfTable["rhoCool"], 
        plusEqOp<scalar>(), 
        rhoCoolant_
    );
    neutroToFluid.mapTgtToSrc(vsfTable["TCool"], plusEqOp<scalar>(), TCoolant_);
    neutroToFluid.mapTgtToSrc
    (
        vsfTable["TStructures"], 
        plusEqOp<scalar>(), 
        TStructures_
    );

    TFuel_.correctBoundaryConditions();
    TCladding_.correctBoundaryConditions();
    rhoCoolant_.correctBoundaryConditions();
    TCoolant_.correctBoundaryConditions();
    TStructures_.correctBoundaryConditions();
    
    phi_ = fvc::flux(neutroToFluid.mapTgtToSrc(vvfTable["U"]));  
}

void Foam::pointKineticNeutronics::getFieldsLiquidFuel
(
    const volVectorFieldTable& vvfTable,
    const volScalarFieldTable& vsfTable, 
    const meshToMesh& neutroToFluid
)
{
    this->getFields
    (
        vvfTable,
        vsfTable, 
        neutroToFluid
    );
}

scalar Foam::pointKineticNeutronics::GEMReactivity()
{
    scalar flow(0);
    forAll(flowFaces_, i)
    {
        flow += phi_[flowFaces_[i]];
    }
    flow = std::abs(flow);

    //- Sync integral flow value across processors
    reduce(flow, sumOp<scalar>()); 

    flowFraction_ = flowFractionRef_*flow/flowRef_;

    //- Specific FFTF relationship between flowFraction and GEM sodium level
    scalar sodiumLevel = 
        (265.0-539504/(2440.13+sqr(flowFraction_*100)))/100 - sodiumLevelRef_;

    Info << "GEMSodiumLevel = " << sodiumLevel << " m" << endl;

    scalar GEMReactivity(0);
    if (GEMReactivityMap_.size() > 1)
    {
        for (int i = 0; i < GEMReactivityMap_.size()-1; i++)
        {
            Pair<scalar> X0(GEMReactivityMap_[i]);
            Pair<scalar> X1(GEMReactivityMap_[i+1]);
            if (sodiumLevel <= X0.first() and sodiumLevel > X1.first())
            {
                scalar m((X1.second()-X0.second())/(X1.first()-X0.first()));
                GEMReactivity = m*(sodiumLevel-X0.first()) + X0.second();
                break;
            }
        }
    }
    
    return GEMReactivity; //flowFraction_;
}

scalar Foam::pointKineticNeutronics::drivelineReactivity()
{
    scalar TDriveline
    (
        fvc::domainIntegrate
        (
            TStructures_*drivelineFeedbackCellField_
        ).value()/
        max
        (
            fvc::domainIntegrate(drivelineFeedbackCellField_).value(),
            SMALL
        )
    );

    //- Absolute drivelineExpansion in meters
    scalar drivelineExpansion
    (
        drivelineExpansionCoeff_*(TDriveline-TDrivelineRef_)
    );
    Info << drivelineExpansionCoeff_ << " " << TDriveline << " " << TDrivelineRef_ << endl;
    Info << "drivelineExpansion = " << drivelineExpansion << " m" << endl;

    //- Interpolate reactivity value from table
    scalar drivelineReactivity(0);
    label N(CRReactivityMap_.size());
    if (N > 1)
    {
        scalar x0(CRReactivityMap_[0].first());
        scalar x1(CRReactivityMap_[N-1].first());

        //- If the map is indexed for ascending parameter values
        if (x0 > x1)
        {
            for (int i = 0; i < N-1; i++)
            {
                Pair<scalar> p0(CRReactivityMap_[i]);
                Pair<scalar> p1(CRReactivityMap_[i+1]);
                if 
                (
                    drivelineExpansion <= p0.first() 
                and drivelineExpansion > p1.first()
                )
                {
                    scalar m((p1.second()-p0.second())/(p1.first()-p0.first()));
                    drivelineReactivity = 
                        m*(drivelineExpansion-p0.first()) + p0.second();
                    break;
                }
            }
        }

        //- If the map is indexed for descending parameter values
        else
        {
            for (int i = 0; i < N-1; i++)
            {
                Pair<scalar> p0(CRReactivityMap_[i]);
                Pair<scalar> p1(CRReactivityMap_[i+1]);
                if 
                (
                    drivelineExpansion > p0.first() 
                and drivelineExpansion <= p1.first()
                )
                {
                    scalar m((p1.second()-p0.second())/(p1.first()-p0.first()));
                    drivelineReactivity = 
                        m*(drivelineExpansion-p0.first()) + p0.second();
                    break;
                }
            }
        }
        
    }

    Info << " ddrivelineReactivity = " << (1e5*drivelineReactivity) << endl;

    return drivelineReactivity;
}

void Foam::pointKineticNeutronics::printInfo() const
{
    typedef HashTable<scalar, word, word::hash> scalarTable;
    
    Info << "    power = " << power_ << " W" << endl;
    Info << "    flowFraction = " << (100*flowFraction_/flowFractionRef_) 
        << " %" << endl;
    Info << "    totalReactivity = " << (1e5*reactivityValues_["total"]) 
        << " pcm" << endl;

    forAllConstIter
    (
        scalarTable,
        reactivityValues_,
        iter
    )
    {
        if (iter.key() != "total")
        {
            scalar value(1e5*iter());
            Info << "    -> " << iter.key() << " = " << value << " pcm" 
                << endl; 
        }
    }
}

// ************************************************************************* //
