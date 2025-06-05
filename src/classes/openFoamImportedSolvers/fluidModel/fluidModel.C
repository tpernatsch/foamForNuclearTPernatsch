/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fluidModel.H"
#include "volFields.H"
#include "fv.H"
#include "fvc.H"

#include "EulerDdtScheme.H"
#include "backwardDdtScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    // defineTypeNameAndDebug(fluidModel, 0);
    // addToRunTimeSelectionTable(physicsModel, fluidModel, physicsModel);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fluidModel::makePisoControl() const
{
    if (!pisoPtr_.empty())
    {
        FatalErrorIn("void Foam::fluidModel::makePisoControl() const")
            << "pointer already set" << abort(FatalError);
    }

    pisoPtr_.set
    (
        new pisoControl
        (
            const_cast<fvMesh&>
            (
                refCast<const fvMesh>(meshRef_)
            )
        )
    );
}


void Foam::fluidModel::makePimpleControl() const
{
    if (!pimplePtr_.empty())
    {
        FatalErrorIn("void Foam::fluidModel::makePimpleControl() const")
            << "pointer already set" << abort(FatalError);
    }

    pimplePtr_.set
    (
        new pimpleControl
        (
            const_cast<fvMesh&>
            (
                refCast<const fvMesh>(meshRef_)
            )
        )
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fluidModel::updateRobinFsiInterface
(
    const volScalarField& p,
    const volVectorField& U,
    surfaceScalarField& phi,
    surfaceScalarField& rAUf
)
{
   
}


void Foam::fluidModel::CourantNo
(
    scalar& CoNum,
    scalar& meanCoNum,
    scalar& velMag
) const
{
    if (meshRef_.nInternalFaces())
    {
        surfaceScalarField magPhi(mag(phi()));

        if (phi().dimensions() == dimVelocity*dimArea*dimDensity)
        {
            Info << "wtf"<<endl;
            const volScalarField& rho =
                meshRef_.lookupObject<volScalarField>("rho");

            magPhi /= fvc::interpolate(rho);
        }

        const surfaceScalarField SfUfbyDelta
        (
            meshRef_.surfaceInterpolation::deltaCoeffs()*magPhi
        );

        const scalar deltaT = mesh_.time().deltaT().value();

        CoNum = max(SfUfbyDelta/meshRef_.magSf()).value()*deltaT;

        meanCoNum = (sum(SfUfbyDelta)/sum(meshRef_.magSf())).value()*deltaT;

        velMag = max(magPhi/meshRef_.magSf()).value();
    }

    Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum
        << " velocity magnitude: " << velMag
        << endl;
}


void Foam::fluidModel::CourantNo() const
{
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
    scalar velMag = 0.0;
    CourantNo(CoNum, meanCoNum, velMag);
}



void Foam::fluidModel::continuityErrs()
{
    const volScalarField contErr(fvc::div(phi()));

    const scalar sumLocalContErr = mesh_.time().deltaT().value()*
        mag(contErr)().weightedAverage(meshRef_.V()).value();

    const scalar globalContErr = mesh_.time().deltaT().value()*
        contErr.weightedAverage(meshRef_.V()).value();

    cumulativeContErr_ += globalContErr;

    Info<< "time step continuity errors : sum local = "
        << sumLocalContErr << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr_
        << endl;
}

void Foam::fluidModel::boundPU
(
    volScalarField& p,
    volVectorField& U
) const
{
    // Bound the pressure
    dimensionedScalar p1 = min(p);
    dimensionedScalar p2 = max(p);

    if (p1 < pMin_ || p2 > pMax_)
    {
        Info<< "p: " << p1.value() << " " << p2.value()
            << ".  Bounding." << endl;

        p.max(pMin_);
        p.min(pMax_);
        p.correctBoundaryConditions();
    }

    // Bound the velocity
    volScalarField magU(mag(U));
    dimensionedScalar U1(max(magU));

    if (U1 > UMax_)
    {
        Info<< "U: " << U1.value() << ".  Bounding." << endl;

        volScalarField Ulimiter
        (
            pos(magU - UMax_)*UMax_/(magU + smallU_) + neg(magU - UMax_)
        );
        Ulimiter.max(scalar(0));
        Ulimiter.min(scalar(1));

        U *= Ulimiter;
        U.correctBoundaryConditions();
    }
}

Foam::meshObjects::gravity Foam::fluidModel::readG() const
{
    // Note: READ_IF_PRESENT is incorreclty implemented within the
    // uniformDimensionedField constructor so we will use a work-around here

    // Check if waveProperties was read from disk; if so, then g must be read
    // from disk
    IOobject wavePropertiesHeader
    (
        "waveProperties",
        mesh_.time().caseConstant(),
        meshRef_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false // do not register
    );

    // Check if g exists on disk
    IOobject gHeader
    (
        "g",
        mesh_.time().caseConstant(),
        meshRef_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false // do not register
    );

    if
    (
        wavePropertiesHeader.typeHeaderOk<IOdictionary>(true)
     && !gHeader.typeHeaderOk<uniformDimensionedVectorField>(true)
    )

    {
        FatalErrorIn("::readG() const")
            << "g field not found in the constant directory: the g field "
            << "must be specified when the waveProperties dictionary is "
            << "present!" << abort(FatalError);
    }

    // The if-else-if structure is broken to keep the compiler happy with the
    // lack of a return statement above

    if
    (
        wavePropertiesHeader.typeHeaderOk<IOdictionary>(true)
     || gHeader.typeHeaderOk<uniformDimensionedVectorField>(true)
    )
    {
        Info<< "Reading g from constant directory" << endl;
        return meshObjects::gravity(mesh_.time());

    }
    else
    {
        Info<< "g field not found in constant directory: initialising to zero"
            << endl;

        return meshObjects::gravity(mesh_.time());

    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidModel::fluidModel
(
    const word& type,
    Time& runTime,
    dynamicFvMesh& mesh,
    const word& region,
    const bool constructNull
)
:
    solver(mesh),
    IOdictionary
    (
        // If region == "region0" then read from the main case
        // Otherwise, read from the region/sub-mesh directory e.g.
        // constant/fluid or constant/solid
        bool(region == dynamicFvMesh::defaultRegion)
      ? IOobject
        (
            "fluidProperties",
            runTime.caseConstant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
      : IOobject
        (
            "fluidProperties",
            runTime.caseConstant(),
            region, // using 'local' property of IOobject
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    meshRef_(mesh),
   /* meshPtr_
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                region,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    ),*/
    fluidProperties_(subDict(type + "Coeffs")),
    pisoPtr_(),
    pimplePtr_(),
    waveProperties_
    (
        IOobject
        (
            "waveProperties",
            runTime.constant(),
            meshRef_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    g_(readG()),
    Uheader_("U", runTime.timeName(), meshRef_, IOobject::MUST_READ),
    pheader_("p", runTime.timeName(), meshRef_, IOobject::MUST_READ),
    UPtr_
    (
        constructNull
      ? nullptr
      : new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                meshRef_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            meshRef_,
            dimensionedVector("zero", dimVelocity, vector::zero)
        )
    ),
    pPtr_
    (
        constructNull
      ? nullptr
      : new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                meshRef_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            meshRef_,
            dimensionedScalar("zero", dimPressure, 0.0)
        )
    ),
    gradUPtr_
    (
        constructNull
      ? nullptr
      : new volTensorField
        (
            IOobject
            (
                "grad(U)",
                runTime.timeName(),
                meshRef_
            ),
            meshRef_,
            dimensionedTensor("zero", dimVelocity/dimLength, tensor::zero)
        )
    ),
    gradpPtr_
    (
        constructNull
      ? nullptr
      : new volVectorField
        (
            IOobject
            (
                "grad(p)",
                runTime.timeName(),
                meshRef_
            ),
            meshRef_,
            dimensionedVector("zero", p().dimensions()/dimLength, vector::zero)
        )
    ),
    phiPtr_
    (
        constructNull
      ? nullptr
      : new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                meshRef_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(U()) & meshRef_.Sf()
        )
    ),
    adjustTimeStep_
    (
        runTime.controlDict().lookupOrDefault<Switch>("adjustTimeStep", false)
    ),
    maxCo_
    (
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0)
    ),
    maxDeltaT_
    (
        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT)
    ),
    pMin_("pMin", dimPressure, 0),
    pMax_("pMax", dimPressure, 0),
    UMax_("UMax", dimVelocity, 0),
    smallU_("smallU", dimVelocity, 1e-10),
    cumulativeContErr_(0.0),
    fvOptions_(fv::options::New(meshRef_)),
    fsiMeshUpdate_(false),
    fsiMeshUpdateChanged_(false)
{
    if (!constructNull)
    {
        gradUPtr_() = fvc::grad(UPtr_());
        gradpPtr_() = fvc::grad(pPtr_());

        pMin_.dimensions().reset(p().dimensions());
        pMax_.dimensions().reset(p().dimensions());
        if (meshRef_.solutionDict().found("fieldBounds"))
        {
            dictionary fieldBounds = meshRef_.solutionDict().subDict("fieldBounds");
            fieldBounds.lookup(p().name())
                >> pMin_.value() >> pMax_.value();
            fieldBounds.lookup(U().name())
                >> UMax_.value();
        }
    }

    // Check if any finite volume option is present
    if (!fvOptions_.optionList::size())
    {
        Info << "No finite volume options present\n" << endl;
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidModel::~fluidModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pisoControl& Foam::fluidModel::piso()
{
    if (pisoPtr_.empty())
    {
        makePisoControl();
    }

    return pisoPtr_();
}


Foam::pimpleControl& Foam::fluidModel::pimple()
{
    if (pimplePtr_.empty())
    {
        makePimpleControl();
    }

    return pimplePtr_();
}


/*Foam::tmp<Foam::vectorField> Foam::fluidModel::faceZoneViscousForce
(
    const label interfaceI
) const
{
    const vectorField patchVF
    (
        patchViscousForce(globalPatches()[interfaceI].patch().index())
    );

    return globalPatches()[interfaceI].patchFaceToGlobal(patchVF);
}


Foam::tmp<Foam::scalarField> Foam::fluidModel::faceZonePressureForce
(
    const label interfaceI
) const
{
    const scalarField patchPF
    (
        patchPressureForce(globalPatches()[interfaceI].patch().index())
    );

    return globalPatches()[interfaceI].patchFaceToGlobal(patchPF);
}


Foam::tmp<Foam::scalarField> Foam::fluidModel::faceZoneTemperature
(
    const label interfaceI
) const
{
    const scalarField patchT
    (
        patchTemperature(globalPatches()[interfaceI].patch().index())
    );

    return globalPatches()[interfaceI].patchFaceToGlobal(patchT);
}


Foam::tmp<Foam::scalarField> Foam::fluidModel::faceZoneHeatFlux
(
    const label interfaceI
) const
{
    const scalarField patchHF
    (
        patchHeatFlux(globalPatches()[interfaceI].patch().index())
    );

    return globalPatches()[interfaceI].patchFaceToGlobal(patchHF);
}


Foam::tmp<Foam::scalarField> Foam::fluidModel::faceZoneHeatTransferCoeff
(
    const label interfaceI
) const
{
    const scalarField patchHTC
    (
        patchHeatTransferCoeff(globalPatches()[interfaceI].patch().index())
    );

    return globalPatches()[interfaceI].patchFaceToGlobal(patchHTC);
}
*/

void Foam::fluidModel::UisRequired()
{
    if (!Uheader_.typeHeaderOk<volVectorField>(true))
    {
        FatalErrorIn("::UisRequired()")
            << "This fluidModel requires the 'U' field to be specified!"
            << abort(FatalError);
    }
}


void Foam::fluidModel::pisRequired()
{
    if (!pheader_.typeHeaderOk<volScalarField>(true))
    {
        FatalErrorIn("::pisRequired()")
            << "This fluidModel requires the 'p' field to be specified!"
            << abort(FatalError);
    }
}



void Foam::fluidModel::setDeltaT(Time& runTime)
{
    if (adjustTimeStep_)
    {
        // Calculate the maximum Courant number
        // Careful to use the relative flux in the calculation
        // We have to be careful when we call makeRelative and makeAbsolute
        scalar CoNum = 0.0;
        scalar meanCoNum = 0.0;
        scalar velMag = 0.0;
        fvc::makeRelative(phi(), U());
        CourantNo(CoNum, meanCoNum, velMag);
        fvc::makeAbsolute(phi(), U());

        scalar maxDeltaTFact = maxCo_/(CoNum + SMALL);
        scalar deltaTFact =
            min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

        runTime.setDeltaT
        (
            min
            (
                deltaTFact*runTime.deltaT().value(),
                maxDeltaT_
            )
        );

        Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
    }
}




// void Foam::fluidModel::writeFields(const Time& runTime)
// {
//     runTime.write();
// }


// bool Foam::fluidModel::read()
// {
//     if (regIOobject::read())
//     {
//         fluidProperties_ = subDict("Coeffs");

//         return true;
//     }
//     else
//     {
//         return false;
//     }
// }

// void Foam::fluidModel::end()
// {
//     this->IOobject::rename(this->IOobject::name()+".withDefaultValues");
//     this->regIOobject::write();
// }
// ************************************************************************* //
